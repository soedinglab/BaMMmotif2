#include "FDR.h"

#include <float.h>		// -FLT_MAX

FDR::FDR( Motif* motif, size_t cvFold ){

	motif_ = motif;
	float fold = ( float )Global::negSequenceSet->getN() / ( float )Global::posSequenceSet->getN();
	mFold_ = ( Global::negSeqGiven ) ? fold : ( float )Global::mFold;
	cvFold_ = cvFold;

	for( size_t k = 0; k < motif_->getK()+8; k++ ){
		Y_.push_back( ipow( Alphabet::getSize(), k ) );
	}

	occ_frac_ = 0.0f;
	occ_mult_ = 0.0f;

}

FDR::~FDR(){

}

void FDR::evaluateMotif( size_t n ){

	std::vector<Sequence*> posSeqs = Global::posSequenceSet->getSequences();
	std::vector<Sequence*> negSeqs = Global::negSequenceSet->getSequences();

	std::vector<std::vector<float>> mops_scores;
	std::vector<float> zoops_scores;
	std::vector<std::vector<float>> scores;

	for( size_t fold = 0; fold < cvFold_; fold++ ){

		if( Global::verbose ){
			fprintf(stderr, " ________________________________\n"
							"|                                |\n"
							"|  Cross validation for fold-%d   |\n"
							"|________________________________|\n\n", (int)fold+1 );
		}

		Motif* motif = new Motif( *motif_ );			// deep copy the initial motif

		// assign training folds
		std::vector<size_t> trainsetFolds;
		trainsetFolds.resize( cvFold_ - 1 );
		for( size_t f = 0; f < cvFold_; f++ ){
			if( f != fold ){
				trainsetFolds.push_back( f );
			}
		}

		// draw sequences for each test set
		std::vector<Sequence*> testSet;
		for( size_t i = 0; i < Global::posFoldIndices[fold].size(); i++ ){
			testSet.push_back( posSeqs[Global::posFoldIndices[fold][i]] );
		}

		// obtain background model for each training set
		BackgroundModel* bgModel;
		if( !Global::bgModelGiven ){
			bgModel = new BackgroundModel( *Global::negSequenceSet,
											Global::bgModelOrder,
											Global::bgModelAlpha,
											Global::interpolateBG,
											Global::negFoldIndices,
											trainsetFolds );

		} else {
			bgModel = new BackgroundModel( Global::bgModelFilename );
		}

		// learn motif from each training set
		if( Global::EM ){
			ModelLearning model( motif, bgModel, trainsetFolds );
			model.EM();
		} else if ( Global::CGS ){
			ModelLearning model( motif, bgModel, trainsetFolds );
			model.GibbsSampling();
		}

		// score positive test sequences with the (learned) motif
		//scores = scoreSequenceSet( motif, bgModel, testSet );
		ScoreSeqSet score_testset( motif, bgModel, testSet );
		score_testset.score();

		if( Global::mops ){
			mops_scores = score_testset.getMopsScores();
			for( size_t n = 0; n < testSet.size(); n++ ){
				posScoreAll_.insert( std::end( posScoreAll_ ), std::begin( mops_scores[n] ), std::end( mops_scores[n] ) );
			}
		}

		if( Global::zoops ){
			zoops_scores = score_testset.getZoopsScores();
			posScoreMax_.insert( std::end( posScoreMax_ ), std::begin( zoops_scores ), std::end( zoops_scores ) );
		}

		if( !Global::negSeqGiven ){
			std::vector<std::unique_ptr<Sequence>> negSet;
			// generate negative sequence set
			SeqGenerator seqs( testSet );
			negSet = seqs.sample_negative_seqset( Global::mFold );
			// score negative sequence set
			scores = scoreSequenceSet( motif, bgModel, negSet );

			if( Global::mops ){
				negScoreAll_.insert( std::end( negScoreAll_ ), std::begin( scores[0] ), std::end( scores[0] ) );
			}
			if( Global::zoops ){
				negScoreMax_.insert( std::end( negScoreMax_ ), std::begin( scores[1] ), std::end( scores[1] ) );
			}
		} else {
			std::vector<Sequence*> negSet;
			for( size_t i = 0; i < Global::negFoldIndices[fold].size(); i++ ){
				negSet.push_back( negSeqs[Global::negFoldIndices[fold][i]] );
			}
			// score negative sequence set
			ScoreSeqSet score_negset( motif, bgModel, negSet );
			score_negset.score();
			if( Global::mops ){
				mops_scores = score_negset.getMopsScores();
				for( size_t n = 0; n < negSet.size(); n++ ){
					negScoreAll_.insert( std::end( negScoreAll_ ), std::begin( mops_scores[n] ), std::end( mops_scores[n] ) );
				}
			}
			if( Global::zoops ){
				zoops_scores = score_negset.getZoopsScores();
				negScoreMax_.insert( std::end( negScoreMax_ ), std::begin( zoops_scores ), std::end( zoops_scores ) );
			}
		}

		if( motif ) 	delete motif;
		if( bgModel ) 	delete bgModel;
	}

	if( Global::savePRs ){
		if( Global::verbose ){
			fprintf(stderr, " __________________________________\n"
							"|                                  |\n"
							"|  calculate precision and recall  |\n"
							"|__________________________________|\n\n" );
		}
		calculatePR();
	}

	if( Global::savePvalues ){
		if( Global::verbose ){
			fprintf(stderr, " ______________________\n"
							"|                      |\n"
							"|  calculate P-values  |\n"
							"|______________________|\n\n" );
		}
		calculatePvalues();
	}

}

void FDR::calculatePR(){

	size_t posN = Global::posSequenceSet->getN();
	size_t negN = ( Global::negSeqGiven ) ? Global::negSequenceSet->getN() : posN * Global::mFold;

	// for MOPS model:
	if( Global::mops ){
		// Sort log odds scores in descending order
		std::sort( posScoreAll_.begin(), posScoreAll_.end(), std::greater<float>() );
		std::sort( negScoreAll_.begin(), negScoreAll_.end(), std::greater<float>() );

		// Rank and score these log odds score values
		size_t idx_posAll = 0;
		size_t idx_negAll = 0;								// index for arrays storing the complete log odds scores
		float E_TP_MOPS = 0.0f;								// maximal distance between TFP and FP = expectation value of TP
		size_t idx_max = posN+negN;							// index when the distance between TFP and FP reaches maximum
															// otherwise, set total number as initial cutoff
		size_t len_all = posScoreAll_.size() + negScoreAll_.size();

		for( size_t i = 0; i < len_all; i++ ){
			if( posScoreAll_[idx_posAll] >= negScoreAll_[idx_negAll] ||
					idx_posAll == 0 || idx_negAll == len_all-1 ){
				idx_posAll++;
			} else {
				idx_negAll++;
			}

			MOPS_TP_.push_back( ( float )idx_posAll - ( float )idx_negAll / mFold_ );
			MOPS_FP_.push_back( ( float )idx_negAll / mFold_ );

			if( E_TP_MOPS == MOPS_TP_[i] ){
				idx_max = i;
	//			break;										// stop when the distance between TFP and FP reaches maximum
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
	if( Global::zoops ){
		PN_Pvalue_.clear();
		// Sort log odds scores in descending order
		std::sort( posScoreMax_.begin(), posScoreMax_.end(), std::greater<float>() );
		std::sort( negScoreMax_.begin(), negScoreMax_.end(), std::greater<float>() );

		// Rank and score these log odds score values
		size_t idx_posMax = 0;
		size_t idx_negMax = 0;
		size_t min_idx_pos = 0;
		size_t posN_est = static_cast<size_t>( Global::q * ( float )posN );

		for( size_t i = 0; i < posN + negN; i++ ){
			if( posScoreMax_[idx_posMax] >= negScoreMax_[idx_negMax] ||
				idx_posMax == 0 || idx_negMax == posN+negN-1 ){
				idx_posMax++;
			} else {
				idx_negMax++;
			}

			ZOOPS_TP_.push_back( ( float )idx_posMax );
			ZOOPS_FP_.push_back( ( float )idx_negMax / mFold_ );
			PN_Pvalue_.push_back( ( ( float )idx_negMax + 0.5f )
					/ ( mFold_ * ( float )posScoreMax_.size() + 1.0f ) );

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

	// Method 1:

	// for MOPS model:
	if( Global::mops ){
		// Sort log odds scores in ascending order
		std::sort( negScoreAll_.begin(), negScoreAll_.end(), std::less<float>() );
		std::sort( posScoreAll_.begin(), posScoreAll_.end(), std::greater<float>() );
		for( size_t i = 0; i < posScoreAll_.size(); i++ ){
			// Return iterator to lower/upper bound
			int low, up;
			low = static_cast<int>( std::distance( negScoreAll_.begin(),
					std::lower_bound( negScoreAll_.begin(), negScoreAll_.end(), posScoreAll_[i] ) ) );
			up =  static_cast<int>( std::distance( negScoreAll_.begin(),
					std::upper_bound( negScoreAll_.begin(), negScoreAll_.end(), posScoreAll_[i] ) ) );
			float p = 1.0f - ( float )( up + low ) / ( 2.0f * ( float ) negScoreAll_.size() );
			// avoid the rounding errors, such as p-value = 0 or p-value > 1
			if( p < 1e-6 ) p = 0.000001f;
			if( p > 1.0f ) p = 1.0f;
			MOPS_Pvalue_.push_back( p );
		}
	}

	// for ZOOPS model:
	if( Global::zoops ){
		// Sort log odds scores in ascending order
		std::sort( negScoreMax_.begin(), negScoreMax_.end(), std::less<float>() );
		std::sort( posScoreMax_.begin(), posScoreMax_.end(), std::greater<float>() );
		for( size_t i = 0; i < posScoreMax_.size(); i++ ){
			// Return iterator to lower/upper bound
			int low, up;
			low = static_cast<int>( std::distance( negScoreMax_.begin(),
					std::lower_bound( negScoreMax_.begin(), negScoreMax_.end(), posScoreMax_[i] ) ) );
			up =  static_cast<int>( std::distance( negScoreMax_.begin(),
					std::upper_bound( negScoreMax_.begin(), negScoreMax_.end(), posScoreMax_[i] ) ) );
			float p = 1.0f - ( float )( up + low ) / ( 2.0f * ( float )negScoreMax_.size() );
			// avoid the rounding errors, such as p-value = 0 or p-value > 1
			if( p < 1e-6 ) p = 0.000001f;
			if( p > 1.0f ) p = 1.0f;
			ZOOPS_Pvalue_.push_back( p );
		}
	}

/*	// Method 2:
	int M;
	if( Global::negSeqGiven ){
		M = Global::negSequenceSet->getN() / Global::posSequenceSet->getN();
	} else {
		M = mFold_;
	}

	// for MOPS model:
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

	// for ZOOPS model:
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
					/ ( ( float )( M * posScoreMax_.size() ) + 1.0f ) );
		} else {
			idx_negMax++;
		}
	}*/

}

// score sequences for the positive sequence set
std::vector<std::vector<float>> FDR::scoreSequenceSet( Motif* motif, BackgroundModel* bg, const std::vector<std::unique_ptr<Sequence>> & seqSet ){

	std::vector<std::vector<float>> scores( 2 );			// scores[0]: store the log odds scores at all positions of each sequence
															// scores[1]: store maximal log odds scores of each sequence
	size_t K = motif->getK();
	size_t W = motif->getW();
	float maxScore;											// maximal logOddsScore over all positions for each sequence

	// pre-calculate log odds scores given motif and bgModel
	motif->calculateS( bg->getV() );
	float** s = motif->getS();

	for( size_t n = 0; n < seqSet.size(); n++ ){

		size_t LW1 = seqSet[n]->getL() - W + 1;
		int* kmer = seqSet[n]->getKmer();

		maxScore = -FLT_MAX;

		for( size_t i = 0; i < LW1; i++ ){

			float logOdds = 0.0f;

			for( size_t j = 0; j < W; j++ ){

				int y = ( kmer[i+j] >= 0 ) ? kmer[i+j] % static_cast<int>( Y_[K+1] ) : -1;

				logOdds += ( y >= 0 ) ? s[y][j] : 0;

			}

			// take all the log odds scores for MOPS model:
			scores[0].push_back( logOdds );

			// take the largest log odds score for ZOOPS model:
			if( logOdds > maxScore ){
				maxScore = ( logOdds > maxScore ) ? logOdds : maxScore;
			}
		}

		scores[1].push_back( maxScore );

	}
	return scores;
}

void FDR::print(){

}

void FDR::write( size_t n ){

	std::string opath = std::string( Global::outputDirectory ) + '/'
			+ Global::posSequenceBasename + "_motif_" + std::to_string( n+1 );

	if( Global::savePRs ){
		/**
		 * save FDR results in two flat files for obtaining AUSFC:
		 * (1) posSequenceBasename.zoops.stats: 	TP, FP, FDR, recall, p-values and mFold for ZOOPS model
		 * (2) posSequenceBasename.mops.stats:		TP, FP, FDR, recall for MOPS model
		 */

		// for ZOOPS model:
		if( Global::zoops ){
			std::string opath_zoops_stats = opath + ".zoops.stats";
			std::ofstream ofile_zoops( opath_zoops_stats );
			// the headers:
			ofile_zoops << "TP" 	<< '\t'
						<< "FP" 	<< '\t'
						<< "FDR" 	<< '\t'
						<< "Recall"	<< '\t'
						<< "p-value"<< '\t'
						<< mFold_ 	<< std::endl;

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
		if( Global::mops ){
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

	if( Global::savePvalues ){
		/**
		 * save FDR results in two flat files for ranking motifs:
		 * (1) posSequenceBasename.zoops.pvalues:	p-values for ZOOPS model
		 * (2) posSequenceBasename.mops.pvalues:	p-values for MOPS model
		 */
		std::string opath = std::string( Global::outputDirectory ) + '/'
				+ Global::posSequenceBasename + "_motif_" + std::to_string( n+1 );

		if( Global::zoops ){
			std::string opath_zoops = opath + ".zoops.pvalues";
			std::ofstream ofile_zoops( opath_zoops );
			for( size_t i = 0; i < ZOOPS_Pvalue_.size(); i++ ){
				ofile_zoops << std::setprecision( 3 )
							<< ZOOPS_Pvalue_[i] << std::endl;
			}
		}

		if( Global::mops ){
			std::string opath_mops = opath + ".mops.pvalues";
			std::ofstream ofile_mops( opath_mops );
			for( size_t i = 0; i < MOPS_Pvalue_.size(); i++ ){
				ofile_mops  << std::setprecision( 3 )
							<< MOPS_Pvalue_[i] << std::endl;
			}
		}
	}

	if( Global::saveLogOdds ){

		/**
		 * save log odds scores into two files for plotting the distribution of log odds scores:
		 * (1) posSequenceBasename.zoops.logOdds
		 * (2) posSequenceBasename.mops.logOdds
		 */

		std::string opath = std::string( Global::outputDirectory ) + '/'
				+ Global::posSequenceBasename + "_motif_" + std::to_string( n+1 );

		if( Global::zoops ){
			std::string opath_zoops_logOdds = opath + ".zoops.logOdds";
			std::ofstream ofile_zoops_logOdds( opath_zoops_logOdds );
			ofile_zoops_logOdds << "positive" << '\t'
								<< "negative" << std::endl;
			for( size_t i = 0; i < posScoreMax_.size(); i++ ){
				ofile_zoops_logOdds	<< std::setprecision( 6 )
									<< posScoreMax_[i] << '\t'
									<< negScoreMax_[i*( size_t )mFold_] << std::endl;
									// * this is not proper for the occasion where negative sequence is given
			}
		}

		if( Global::mops ){
			std::string opath_mops_logOdds = opath + ".mops.logOdds";
			std::ofstream ofile_mops_logOdds( opath_mops_logOdds );
			ofile_mops_logOdds 	<< "positive" << '\t'
								<< "negative" << std::endl;
			for( size_t i = 0; i < posScoreAll_.size(); i++ ){
				ofile_mops_logOdds 	<< std::setprecision( 6 )
									<< posScoreAll_[i] << '\t'
									<< negScoreAll_[i*( size_t )mFold_] << std::endl;
									// * this is not proper for the occasion where negative sequence is given
			}
		}
	}
}
