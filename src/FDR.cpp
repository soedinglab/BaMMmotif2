#include "FDR.h"

FDR::FDR( std::vector<Sequence*> posSeqs,
			std::vector<Sequence*> negSeqs,
			float q,
			Motif* motif,
			size_t cvFold ){

	posSeqs_	= posSeqs;
	negSeqs_	= negSeqs;
	q_ 			= q;
	motif_ 		= motif;
	cvFold_		= cvFold;
	occ_frac_	= 0.0f;
	occ_mult_	= 0.0f;

}

FDR::~FDR(){

}

void FDR::evaluateMotif(){

	std::vector<std::vector<float>> mops_scores;
	std::vector<float> 				zoops_scores;

	// generate sequences from positive sequences with masked motif
	// deep copy the initial motif
	Motif* motif_opti = new Motif( *motif_ );
	BackgroundModel* bgModel_opti = new BackgroundModel( posSeqs_,
									Global::bgModelOrder,
									Global::bgModelAlpha );
	ModelLearning model_opti( motif_opti, bgModel_opti, posSeqs_, q_ );
	model_opti.EM();
	float** r = model_opti.getR();
	SeqGenerator artificalset( posSeqs_, motif_opti );
	std::vector<std::unique_ptr<Sequence>> B1SeqSetPrime;
	B1SeqSetPrime = artificalset.arti_negset_motif_masked( r );
	delete motif_opti;
	delete bgModel_opti;

	/**
	 * Cross Validation
	 */
	for( size_t fold = 0; fold < cvFold_; fold++ ){

		// deep copy the initial motif
		Motif* motif = new Motif( *motif_ );

		/**
		 * Draw sequences for each training and test set
		 */
		std::vector<Sequence*> testSet;
		std::vector<Sequence*> trainSet;
		for( size_t n = 0; n < posSeqs_.size()-cvFold_; n+=cvFold_ ){
			for( size_t f = 0; f < cvFold_; f++ ){
				if( f != fold ){
					trainSet.push_back( posSeqs_[n+f] );
				} else {
					testSet.push_back( posSeqs_[n+f] );
				}
			}
		}

		/**
		 * Generate negative sequence set
		 * for generating bg model and scoring
		 */

		std::vector<Sequence*> B1set;
		std::vector<Sequence*> B2set;
		std::vector<Sequence*> B3set;
		std::vector<Sequence*> B1setPrime;
		// draw B3set from the given negative sequence set
		for( size_t n = 0; n < negSeqs_.size()-cvFold_; n+=cvFold_ ){
			B3set.push_back( negSeqs_[n+fold] );
		}

		// sample negative sequence set B1set based on s-mer frequencies
		// from positive training sequence set
		std::vector<std::unique_ptr<Sequence>> B1Seqs;
		SeqGenerator b1seqs( trainSet );
		B1Seqs = b1seqs.arti_negset( Global::mFold );
		// convert unique_ptr to regular pointer
		for( size_t n = 0; n < B1Seqs.size(); n++ ){
			B1set.push_back( B1Seqs[n].release() );
		}

		// sample negative sequence set B2set based on s-mer frequencies
		// from the given sampled negative sequence set
		std::vector<std::unique_ptr<Sequence>> B2Seqs;
		SeqGenerator b2seqs( B3set );
		B2Seqs = b2seqs.arti_negset( 1 );
		// convert unique_ptr to regular pointer
		for( size_t n = 0; n < B2Seqs.size(); n++ ){
			B2set.push_back( B2Seqs[n].release() );
		}

		// draw B1setPrime from the positive sequence set with masked motifs
		for( size_t n = 0; n < B1SeqSetPrime.size()-cvFold_; n+=cvFold_ ){
			B1setPrime.push_back( B1SeqSetPrime[n+fold].release() );
		}

		/**
		 * Generate background model from negative set
		 */
		BackgroundModel* bgModel;

		if( Global::bgModelGiven ){
			// get the given bg model
			bgModel = new BackgroundModel( Global::bgModelFilename );
		} else {
			// generate K-th background model from 20% of negative sequences
			bgModel = new BackgroundModel( B1set,
											Global::bgModelOrder,
											Global::bgModelAlpha );
		}

		/**
		 * Training
		 */
		// learn motif from each training set
		if( EM_ ){
			ModelLearning model( motif, bgModel, trainSet, q_ );
			model.EM();
		} else if ( CGS_ ){
			ModelLearning model( motif, bgModel, trainSet, q_ );
			model.GibbsSampling();
		}

		/**
		 * Testing
		 */
		// score positive test sequences with (learned) motif
		ScoreSeqSet score_testset( motif, bgModel, testSet );
		score_testset.score();

		if( mops_ ){
			mops_scores = score_testset.getMopsScores();
			for( size_t n = 0; n < testSet.size(); n++ ){
				posScoreAll_.insert( std::end( posScoreAll_ ),
						std::begin( mops_scores[n] ),
						std::end( mops_scores[n] ) );
			}
		}

		if( zoops_ ){
			zoops_scores = score_testset.getZoopsScores();
			posScoreMax_.insert( std::end( posScoreMax_ ),
					std::begin( zoops_scores ),
					std::end( zoops_scores ) );
		}

		// score negative sequence set
		ScoreSeqSet score_negset( motif, bgModel, B3set );
		score_negset.score();
		if( mops_ ){
			mops_scores = score_negset.getMopsScores();
			for( size_t n = 0; n < B3set.size(); n++ ){
				negScoreAll_.insert( std::end( negScoreAll_ ),
						std::begin( mops_scores[n] ),
						std::end( mops_scores[n] ) );
			}
		}
		if( zoops_ ){
			zoops_scores = score_negset.getZoopsScores();
			negScoreMax_.insert( std::end( negScoreMax_ ),
					std::begin( zoops_scores ),
					std::end( zoops_scores ) );
		}

		if( motif ) 				delete motif;
		if( bgModel )				delete bgModel;


	}

	if( savePRs_ ){

		fprintf(stderr, " __________________________________\n"
						"|                                  |\n"
						"|  calculate precision and recall  |\n"
						"|__________________________________|\n\n" );

		calculatePR();
	}

	if( savePvalues_ ){

		fprintf(stderr, " ______________________\n"
						"|                      |\n"
						"|  calculate P-values  |\n"
						"|______________________|\n\n" );

		calculatePvalues();
	}

}

void FDR::calculatePR(){

	size_t posN = posSeqs_.size();
	size_t negN = negScoreMax_.size();
	float mFold = ( float )negN / ( float )posN;
	// for MOPS model:
	if( mops_ ){
		// Sort log odds scores in descending order
		std::sort( posScoreAll_.begin(), posScoreAll_.end(), std::greater<float>() );
		std::sort( negScoreAll_.begin(), negScoreAll_.end(), std::greater<float>() );

		// Rank and score these log odds score values
		size_t idx_posAll = 0;
		size_t idx_negAll = 0;		// index for arrays storing log odds scores
		float E_TP_MOPS = 0.0f;		// maximal distance between TFP and
									// FP = expectation value of TP
		size_t idx_max = posN+negN;	// index when the distance between TFP and
									// FP reaches maximum; otherwise, set the
									// the number as initial cutoff
		size_t len_all = posScoreAll_.size() + negScoreAll_.size();

		for( size_t i = 0; i < len_all; i++ ){
			if( posScoreAll_[idx_posAll] >= negScoreAll_[idx_negAll] ||
					idx_posAll == 0 || idx_negAll == len_all ){
				idx_posAll++;
			} else {
				idx_negAll++;
			}

			MOPS_TP_.push_back( ( float )idx_posAll - ( float )idx_negAll / mFold );
			MOPS_FP_.push_back( ( float )idx_negAll / mFold );

			if( E_TP_MOPS == MOPS_TP_[i] ){
				idx_max = i;
	//			break;				// stop when the distance between TFP and
									// FP reaches maximum
			}
			if( E_TP_MOPS < MOPS_TP_[i] ){
				E_TP_MOPS = MOPS_TP_[i];
			}
		}

		for( size_t i = 0; i < idx_max; i++ ){
			MOPS_FDR_.push_back( MOPS_FP_[i] / ( MOPS_TP_[i] + MOPS_FP_[i] ) );
			MOPS_Rec_.push_back( MOPS_TP_[i] / E_TP_MOPS );
		}

		// the number of motif occurrences per sequence
		occ_mult_ = E_TP_MOPS / ( float )posN;
	}

	// for ZOOPS model:
	if( zoops_ ){
		PN_Pvalue_.clear();
		// Sort log odds scores in descending order
		std::sort( posScoreMax_.begin(), posScoreMax_.end(), std::greater<float>() );
		std::sort( negScoreMax_.begin(), negScoreMax_.end(), std::greater<float>() );

		// Rank and score these log odds score values
		size_t idx_posMax = 0;
		size_t idx_negMax = 0;
		size_t min_idx_pos = 0;
		size_t posN_est = static_cast<size_t>( q_ * ( float )posN );

		for( size_t i = 0; i < posN + negN; i++ ){
			if( posScoreMax_[idx_posMax] >= negScoreMax_[idx_negMax] ||
				idx_posMax == 0 || idx_negMax == posN+negN-1 ){
				idx_posMax++;
			} else {
				idx_negMax++;
			}

			ZOOPS_TP_.push_back( ( float )idx_posMax );
			ZOOPS_FP_.push_back( ( float )idx_negMax / mFold );
			PN_Pvalue_.push_back( ( ( float )idx_negMax + 0.5f )
					/ ( mFold * ( float )posScoreMax_.size() + 1.0f ) );

			// take the faction of q sequences as real positives
			if( idx_posMax == posN_est ){
				min_idx_pos = i + 1;
			}

			ZOOPS_FDR_.push_back( ZOOPS_FP_[i] / ( ZOOPS_TP_[i] + ZOOPS_FP_[i] ) );
			// ZOOPS_FDR_.push_back( 1.0f - ( float )idx_posMax / ( float )(i+1) );
			// ZOOPS_FDR_.push_back( ZOOPS_FP_[i] / ( float )idx_posMax );
			ZOOPS_Rec_.push_back( ZOOPS_TP_[i] / ( float )posN );

		}

		// the fraction of motif occurrence
		occ_frac_ = 1 - ZOOPS_FP_[min_idx_pos] / ( float )posN;
	}
}

void FDR::calculatePvalues(){

	/**
	 *  calculate p-values for the log odds scores from positive sequence set
	 * 	later used for ranking models using fdrtool R package
	 */

	bool use_method_1 = true;
	if( use_method_1 ){
	// Method 1:
	// for MOPS model:
	if( mops_ ){
		// Sort log odds scores in ascending order
		std::sort( negScoreAll_.begin(), negScoreAll_.end(), std::less<float>() );
		std::sort( posScoreAll_.begin(), posScoreAll_.end(), std::less<float>() );
		for( size_t i = 0; i < posScoreAll_.size(); i++ ){
			// Return iterator to lower/upper bound
			int low, up;
			low = static_cast<int>( std::distance( negScoreAll_.begin(),
										std::lower_bound( negScoreAll_.begin(),
												negScoreAll_.end(),
												posScoreAll_[i] ) ) );
			up =  static_cast<int>( std::distance( negScoreAll_.begin(),
										std::upper_bound( negScoreAll_.begin(),
												negScoreAll_.end(),
												posScoreAll_[i] ) ) );
			float p = 1.0f - ( float )( up + low ) /
					( 2.0f * ( float ) negScoreAll_.size() );
			// avoid the rounding errors, such as p-value = 0 or p-value > 1
			if( p < 1e-6 ) p = 0.000001f;
			if( p > 1.0f ) p = 1.0f;
			MOPS_Pvalue_.push_back( p );
		}
	}

	// for ZOOPS model:
	if( zoops_ ){
		// Sort log odds scores in ascending order
		std::sort( negScoreMax_.begin(), negScoreMax_.end(), std::less<float>() );
		std::sort( posScoreMax_.begin(), posScoreMax_.end(), std::less<float>() );
		for( size_t i = 0; i < posScoreMax_.size(); i++ ){
			// Return iterator to lower/upper bound
			int low, up;
			low = static_cast<int>( std::distance( negScoreMax_.begin(),
										std::lower_bound( negScoreMax_.begin(),
												negScoreMax_.end(),
												posScoreMax_[i] ) ) );
			up =  static_cast<int>( std::distance( negScoreMax_.begin(),
										std::upper_bound( negScoreMax_.begin(),
												negScoreMax_.end(),
												posScoreMax_[i] ) ) );
			float p = 1.0f - ( float )( up + low )
					/ ( 2.0f * ( float )negScoreMax_.size() );
			// avoid the rounding errors, such as p-value = 0 or p-value > 1
			if( p < 1e-6 ) p = 0.000001f;
			if( p > 1.0f ) p = 1.0f;
			ZOOPS_Pvalue_.push_back( p );
		}
	}

	} else {
	// Method 2:
	// for MOPS model:
	if( mops_ ){
		// Sort log odds scores in descending order
		std::sort( posScoreAll_.begin(), posScoreAll_.end(), std::greater<float>() );
		std::sort( negScoreAll_.begin(), negScoreAll_.end(), std::greater<float>() );

		// Rank and score these log odds score values
		size_t idx_posAll = 0;
		size_t idx_negAll = 0;
		size_t len_all = posScoreAll_.size() + negScoreAll_.size();

		for( size_t i = 0; i < len_all; i++ ){
			if( posScoreAll_[idx_posAll] >= negScoreAll_[idx_negAll]
				|| idx_posAll == 0 || idx_negAll == len_all - 1 ){
				idx_posAll++;
				MOPS_Pvalue_.push_back( ( ( float )idx_negAll + 0.5f )
						/ ( ( float ) negScoreAll_.size()  + 1.0f ) );
			} else {
				idx_negAll++;
			}

		}
	}

	// for ZOOPS model:
	if( zoops_ ){
		// Sort log odds scores in descending order
		std::sort( posScoreMax_.begin(), posScoreMax_.end(), std::greater<float>() );
		std::sort( negScoreMax_.begin(), negScoreMax_.end(), std::greater<float>() );

		// Rank and score these log odds score values
		size_t idx_posMax = 0;
		size_t idx_negMax = 0;
		size_t len_max = posScoreMax_.size() + negScoreMax_.size();
		for( size_t i = 0; i < len_max; i++ ){
			if( posScoreMax_[idx_posMax] >= negScoreMax_[idx_negMax]
				|| idx_posMax == 0  || idx_negMax == len_max - 1 ){
				idx_posMax++;
				ZOOPS_Pvalue_.push_back( ( ( float )idx_negMax + 0.5f )
						/ ( ( float )negSeqs_.size() + 1.0f ) );
			} else {
				idx_negMax++;
			}
		}
	}

	}
}

void FDR::print(){

}

void FDR::write( char* odir, std::string basename, size_t n ){

	std::string opath = std::string( odir ) + '/' + basename +
						"_motif_" + std::to_string( n );

	if( savePRs_ ){
		/**
		 * save FDR results in two flat files for obtaining AUSFC:
		 * (1) posSequenceBasename.zoops.stats:
		 * TP, FP, FDR, recall, p-values and mFold for ZOOPS model
		 * (2) posSequenceBasename.mops.stats:
		 * TP, FP, FDR, recall for MOPS model
		 */

		// for ZOOPS model:
		if( zoops_ ){
			std::string opath_zoops_stats = opath + ".zoops.stats";
			std::ofstream ofile_zoops( opath_zoops_stats );
			// the headers:
			ofile_zoops << "TP" 	<< '\t'
						<< "FP" 	<< '\t'
						<< "FDR" 	<< '\t'
						<< "Recall"	<< '\t'
						<< "p-value"<< '\t'
						<< ( float )negSeqs_.size() / ( float )posSeqs_.size()
						<< std::endl;

			for( size_t i = 0; i < ZOOPS_FDR_.size(); i++ ){
				ofile_zoops << ZOOPS_TP_[i]  << '\t'
							<< ZOOPS_FP_[i]  << '\t'
							<< ZOOPS_FDR_[i] << '\t'
							<< ZOOPS_Rec_[i] << '\t'
							<< PN_Pvalue_[i] << '\t'
							<< std::endl;
			}
		}

		// for MOPS model:
		if( mops_ ){
			std::string opath_mops_stats = opath + ".mops.stats";
			std::ofstream ofile_mops( opath_mops_stats );
			ofile_mops  << "TP" 	<< '\t'
						<< "FP" 	<< '\t'
						<< "FDR" 	<< '\t'
						<< "Recall"	<< std::endl;

			for( size_t i = 0; i < MOPS_FDR_.size(); i++ ){
				ofile_mops  << MOPS_TP_[i]  << '\t'
							<< MOPS_FP_[i]  << '\t'
							<< MOPS_FDR_[i] << '\t'
							<< MOPS_Rec_[i] << '\t' << std::endl;
			}
		}
	}

	if( savePvalues_ ){
		/**
		 * save FDR results in two flat files for ranking motifs:
		 * (1) posSequenceBasename.zoops.pvalues:	p-values for ZOOPS model
		 * (2) posSequenceBasename.mops.pvalues:	p-values for MOPS model
		 */
		if( zoops_ ){
			std::string opath_zoops = opath + ".zoops.pvalues";
			std::ofstream ofile_zoops( opath_zoops );
			for( size_t i = 0; i < ZOOPS_Pvalue_.size(); i++ ){
				ofile_zoops << std::setprecision( 3 )
							<< ZOOPS_Pvalue_[i] << std::endl;
			}
		}

		if( mops_ ){
			std::string opath_mops = opath + ".mops.pvalues";
			std::ofstream ofile_mops( opath_mops );
			for( size_t i = 0; i < MOPS_Pvalue_.size(); i++ ){
				ofile_mops  << std::setprecision( 3 )
							<< MOPS_Pvalue_[i] << std::endl;
			}
		}
	}

	if( saveLogOdds_ ){

		/**
		 * save log odds scores into two files for
		 * plotting the distribution of log odds scores:
		 * (1) posSequenceBasename.zoops.logOdds
		 * (2) posSequenceBasename.mops.logOdds
		 */
		if( zoops_ ){
			std::string opath_zoops_logOdds = opath + ".zoops.logOdds";
			std::ofstream ofile_zoops_logOdds( opath_zoops_logOdds );
			ofile_zoops_logOdds << "positive" << '\t'
								<< "negative" << std::endl;
			for( size_t i = 0; i < posScoreMax_.size(); i++ ){
				ofile_zoops_logOdds	<< std::setprecision( 6 )
									<< posScoreMax_[i] << '\t'
									<< negScoreMax_[i*negSeqs_.size()/posSeqs_.size()]
									<< std::endl;
			}
		}

		if( mops_ ){
			std::string opath_mops_logOdds = opath + ".mops.logOdds";
			std::ofstream ofile_mops_logOdds( opath_mops_logOdds );
			ofile_mops_logOdds 	<< "positive" << '\t'
								<< "negative" << std::endl;
			for( size_t i = 0; i < posScoreAll_.size(); i++ ){
				ofile_mops_logOdds 	<< std::setprecision( 6 )
									<< posScoreAll_[i] << '\t'
									<< negScoreAll_[i*negSeqs_.size()/posSeqs_.size()]
									<< std::endl;
			}
		}
	}
}
