#include "FDR.h"

#include <float.h>		// -FLT_MAX

FDR::FDR( Motif* motif ){

	motif_ = motif;

	for( int k = 0; k < Global::Yk; k++ ){
		Y_.push_back( ipow( Alphabet::getSize(), k ) );
	}

	occ_frac_ = 0.0f;

	occ_mult_ = 0.0f;

}

FDR::~FDR(){

}

void FDR::evaluateMotif( int n ){

	std::vector<Sequence*> posSeqs = Global::posSequenceSet->getSequences();
	std::vector<Sequence*> negSeqs = Global::negSequenceSet->getSequences();

	std::vector<std::vector<float>> scores;

	int num = Global::posSequenceSet->getN();
	if( num < 5000 ){
		Global::mFold = 50000 / num;
	}

	for( int fold = 0; fold < Global::cvFold; fold++ ){

		if( Global::verbose ){
			fprintf(stderr, " ________________________________\n"
							"|                                |\n"
							"|  Cross validation for fold-%d   |\n"
							"|________________________________|\n\n", fold+1 );
		}

		Motif* motif = new Motif( *motif_ );			// deep copy the initial motif

		// assign training folds
		std::vector<int> trainsetFolds;
		trainsetFolds.resize( Global::cvFold - 1 );
		for( int f = 0; f < Global::cvFold; f++ ){
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
		BackgroundModel* bgModel = new BackgroundModel( *Global::negSequenceSet,
													Global::bgModelOrder,
													Global::bgModelAlpha,
													Global::interpolateBG,
													Global::negFoldIndices,
													trainsetFolds );

		// learn motif from each training set
		if( Global::EM ){
			ModelLearning model( motif, bgModel, trainsetFolds );
			model.EM();
		} else if ( Global::CGS ){
			ModelLearning model( motif, bgModel, trainsetFolds );
			model.GibbsSampling();
		}

		// score positive test sequences with the (learned) motif
		scores = scoreSequenceSet( motif, bgModel, testSet );
		posScoreAll_.insert( std::end( posScoreAll_ ), std::begin( scores[0] ), std::end( scores[0] ) );
		posScoreMax_.insert( std::end( posScoreMax_ ), std::begin( scores[1] ), std::end( scores[1] ) );

		if( !Global::negSeqGiven ){
			std::vector<std::unique_ptr<Sequence>> negSet;
			// generate negative sequence set
			SeqGenerator seqs( testSet );
			negSet = seqs.sample_negative_seqset();
			// score negative sequence set
			scores = scoreSequenceSet( motif, bgModel, negSet );
		} else {
			std::vector<Sequence*> negSet;
			for( size_t i = 0; i < Global::negFoldIndices[fold].size(); i++ ){
				negSet.push_back( negSeqs[Global::negFoldIndices[fold][i]] );
			}
			// score negative sequence set
			scores = scoreSequenceSet( motif, bgModel, negSet );
		}

		negScoreAll_.insert( std::end( negScoreAll_ ), std::begin( scores[0] ), std::end( scores[0] ) );
		negScoreMax_.insert( std::end( negScoreMax_ ), std::begin( scores[1] ), std::end( scores[1] ) );

		delete motif;
		delete bgModel;
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

// score sequences for the positive sequence set
std::vector<std::vector<float>> FDR::scoreSequenceSet( Motif* motif, BackgroundModel* bg, const std::vector<std::unique_ptr<Sequence>> & seqSet ){

	std::vector<std::vector<float>> scores( 2 );			// scores[0]: store the log odds scores at all positions of each sequence
															// scores[1]: store maximal log odds scores of each sequence
	int K = Global::modelOrder;
	int K_bg = Global::bgModelOrder;
	int W = motif->getW();
	float maxScore;											// maximal logOddsScore over all positions for each sequence

	for( size_t n = 0; n < seqSet.size(); n++ ){

		int LW1 = seqSet[n]->getL() - W + 1;
		int* kmer = seqSet[n]->getKmer();

		maxScore = -FLT_MAX;

		std::vector<float> logOdds( LW1 );					// calculate log odds scores for all the positions

		for( int i = 0; i < LW1; i++ ){

			logOdds[i] = 0.0f;

			for( int j = 0; j < W; j++ ){

				int y = kmer[i+j] % Y_[(i+j < K ) ? i+j+1 : K+1];

				int y_bg = y % Y_[K_bg+1];

				if( y >= 0 ){
					logOdds[i] += ( logf( motif->getV()[K][y][j] + 0.000001f ) - logf( bg->getV()[std::min( K, K_bg )][y_bg] ) );
				}
			}

			// take all the log odds scores for MOPS model:
			scores[0].push_back( logOdds[i] );

			// take the largest log odds score for ZOOPS model:
			if( logOdds[i] > maxScore ){
				maxScore = logOdds[i];
			}
		}

		scores[1].push_back( maxScore );

	}
	return scores;
}

// score sequences for the negative sequence set
std::vector<std::vector<float>> FDR::scoreSequenceSet( Motif* motif, BackgroundModel* bg, std::vector<Sequence*> seqSet ){

	std::vector<std::vector<float>> scores( 2 );			// scores[0]: store the log odds scores at all positions of each sequence
															// scores[1]: store maximal log odds scores of each sequence
	int K = Global::modelOrder;
	int K_bg = Global::bgModelOrder;
	int W = motif->getW();
	float maxScore;											// maximal logOddsScore over all positions for each sequence

	for( size_t n = 0; n < seqSet.size(); n++ ){

		int LW1 = seqSet[n]->getL() - W + 1;
		int* kmer = seqSet[n]->getKmer();

		maxScore = -FLT_MAX;

		std::vector<float> logOdds( LW1 );					// calculate log odds scores for all the positions

		for( int i = 0; i < LW1; i++ ){

			logOdds[i] = 0.0f;

			for( int j = 0; j < W; j++ ){
				int y = kmer[i+j] % Y_[(i+j < K ) ? i+j+1 : K+1];
				int y_bg = y % Y_[K_bg+1];

				if( y >= 0 ){
					logOdds[i] += ( logf( motif->getV()[K][y][j] + 0.000001f ) - logf( bg->getV()[std::min( K, K_bg )][y_bg] ) );
				}
			}

			// take all the log odds scores for MOPS model:
			scores[0].push_back( logOdds[i] );

			// take the largest log odds score for ZOOPS model:
			if( logOdds[i] > maxScore ){
				maxScore = logOdds[i];
			}
		}

		scores[1].push_back( maxScore );

	}
	return scores;
}

void FDR::calculatePR(){

	int M;
	int posN = Global::posSequenceSet->getN();
	int negN;
	if( Global::negSeqGiven ){
		negN = Global::negSequenceSet->getN();
		M = negN / posN;
	} else {
		M = Global::mFold;
		negN = posN * M;

	}

	PN_Pvalue_.clear();

	// for MOPS model:
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

		MOPS_TP_.push_back( ( float )idx_posAll - ( float )idx_negAll / ( float )M );
		MOPS_FP_.push_back( ( float )idx_negAll / ( float )M );

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

	// for ZOOPS model:
	// Sort log odds scores in descending order
	std::sort( posScoreMax_.begin(), posScoreMax_.end(), std::greater<float>() );
	std::sort( negScoreMax_.begin(), negScoreMax_.end(), std::greater<float>() );

	// Rank and score these log odds score values
	int idx_posMax = 0;
	int idx_negMax = 0;
	int min_idx_pos = 0;
	int posN_est = static_cast<int>( Global::q * ( float )posN );

	for( int i = 0; i < posN + negN; i++ ){
		if( posScoreMax_[idx_posMax] >= negScoreMax_[idx_negMax] ||
			idx_posMax == 0 ||  idx_negMax == posN+negN-1 ){
			idx_posMax++;
		} else {
			idx_negMax++;
		}

		ZOOPS_TP_.push_back( ( float )idx_posMax );
		ZOOPS_FP_.push_back( ( float )idx_negMax / ( float )M );
		PN_Pvalue_.push_back( ( ( float )idx_negMax + 0.5f )
				/ ( ( float )( M * posScoreMax_.size() ) + 1.0f ) );

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

void FDR::calculatePvalues(){

	/**
	 *  calculate p-values for the log odds scores from positive sequence set
	 * 	later used for evaluating models using fdrtool R package
	 */

	MOPS_Pvalue_.clear();
	ZOOPS_Pvalue_.clear();

	// Method 1:
	// for MOPS model:
	// Sort log odds scores in ascending order
	std::sort( negScoreAll_.begin(), negScoreAll_.end(), std::less<float>() );
	std::sort( posScoreAll_.begin(), posScoreAll_.end(), std::greater<float>() );
	for( size_t i = 0; i < posScoreAll_.size(); i++ ){
		// Return iterator to lower/upper bound
		size_t low, up;
		low = std::distance( negScoreAll_.begin(), std::lower_bound( negScoreAll_.begin(), negScoreAll_.end(), posScoreAll_[i] ) );
		up =  std::distance( negScoreAll_.begin(), std::upper_bound( negScoreAll_.begin(), negScoreAll_.end(), posScoreAll_[i] ) );
		float p = 1.0f - ( float )( up + low ) / ( 2.0f * ( float ) negScoreAll_.size() );
		// avoid the rounding errors, such as p-value = 0 or p-value > 1
		if( p < 1e-6 ) p = 0.000001f;
		if( p > 1.0f ) p = 1.0f;
		MOPS_Pvalue_.push_back( p );
	}

	// for ZOOPS model:
	// Sort log odds scores in ascending order
	std::sort( negScoreMax_.begin(), negScoreMax_.end(), std::less<float>() );
	std::sort( posScoreMax_.begin(), posScoreMax_.end(), std::greater<float>() );
	for( size_t i = 0; i < posScoreMax_.size(); i++ ){
		// Return iterator to lower/upper bound
		size_t low, up;
		low = std::distance( negScoreMax_.begin(), std::lower_bound( negScoreMax_.begin(), negScoreMax_.end(), posScoreMax_[i] ) );
		up =  std::distance( negScoreMax_.begin(), std::upper_bound( negScoreMax_.begin(), negScoreMax_.end(), posScoreMax_[i] ) );
		float p = 1.0f - ( float )( up + low ) / ( 2.0f * ( float )negScoreMax_.size() );
		// avoid the rounding errors, such as p-value = 0 or p-value > 1
		if( p < 1e-6 ) p = 0.000001f;
		if( p > 1.0f ) p = 1.0f;
		ZOOPS_Pvalue_.push_back( p );
	}

/*	// Method 2:
	int M;
	if( Global::negSeqGiven ){
		M = Global::negSequenceSet->getN() / Global::posSequenceSet->getN();
	} else {
		M = Global::mFold;
	}

	// for MOPS model:
	// Sort log odds scores from large to small
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
	// Sort log odds scores from large to small
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

void FDR::print(){

}

void FDR::write( int n ){

	std::string opath = std::string( Global::outputDirectory ) + '/'
			+ Global::posSequenceBasename + "_motif_" + std::to_string( n+1 );

	if( Global::savePRs ){
		/**
		 * save FDR results in two flat files for obtaining AUSFC:
		 * (1) posSequenceBasename.zoops.stats: 	TP, FP, FDR, recall, p-values and mFold for ZOOPS model
		 * (2) posSequenceBasename.mops.stats:		TP, FP, FDR, recall for MOPS model
		 */
		std::string opath_zoops_stats = opath + ".zoops.stats";
//		std::string opath_mops_stats = opath + ".mops.stats";

		std::ofstream ofile_zoops( opath_zoops_stats );
//		std::ofstream ofile_mops( opath_mops_stats );
		size_t i;

		// for ZOOPS model:
		// the headers:
		ofile_zoops << "TP" << '\t'
					<< "FP" << '\t'
					<< "FDR" << '\t'
					<< "Recall" << '\t'
					<< "p-value" << '\t'
					<< Global::mFold << std::endl;

		for( i = 0; i < ZOOPS_FDR_.size(); i++ ){
			ofile_zoops << ZOOPS_TP_[i]  << '\t'
						<< ZOOPS_FP_[i]  << '\t'
						<< ZOOPS_FDR_[i] << '\t'
						<< ZOOPS_Rec_[i] << '\t'
						<< PN_Pvalue_[i] << std::endl;
		}

/*		// for MOPS model:
		ofile_mops  << "TP" << '\t'
					<< "FP" << '\t'
					<< "FDR" << '\t'
					<< "Recall" << std::endl;


		for( i = 0; i < MOPS_FDR_.size(); i++ ){
			ofile_mops  << MOPS_TP_[i]  << '\t'
						<< MOPS_FP_[i]  << '\t'
						<< MOPS_FDR_[i] << '\t'
						<< MOPS_Rec_[i] << '\t' << std::endl;
		}*/
	}

	if( Global::savePvalues ){
		/**
		 * save FDR results in two flat files for ranking motifs:
		 * (1) posSequenceBasename.zoops.pvalues:	p-values for ZOOPS model
		 * (2) posSequenceBasename.mops.pvalues:	p-values for MOPS model
		 */
		std::string opath = std::string( Global::outputDirectory ) + '/'
				+ Global::posSequenceBasename + "_motif_" + std::to_string( n+1 );

		std::string opath_zoops = opath + ".zoops.pvalues";
		std::string opath_mops = opath + ".mops.pvalues";

		std::ofstream ofile_zoops( opath_zoops );
		std::ofstream ofile_mops( opath_mops );

		for( size_t i = 0; i < ZOOPS_Pvalue_.size(); i++ ){
			ofile_zoops << std::setprecision(3)
						<< ZOOPS_Pvalue_[i] << std::endl;
		}

		for( size_t i = 0; i < MOPS_Pvalue_.size(); i++ ){
			ofile_mops  << std::setprecision(3)
						<< MOPS_Pvalue_[i] << std::endl;
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

		std::string opath_zoops_logOdds = opath + ".zoops.logOdds";
		std::string opath_mops_logOdds = opath + ".mops.logOdds";
		std::ofstream ofile_zoops_logOdds( opath_zoops_logOdds );
		std::ofstream ofile_mops_logOdds( opath_mops_logOdds );

		size_t i;
		int M;
		if( Global::negSeqGiven ){
			M = Global::negSequenceSet->getN() / Global::posSequenceSet->getN();
		} else {
			M = Global::mFold;
		}

		ofile_zoops_logOdds << "positive" << '\t'
							<< "negative" << std::endl;
		for( i = 0; i < posScoreMax_.size(); i++ ){
			ofile_zoops_logOdds	<< std::setprecision(6)
								<< posScoreMax_[i] << '\t'
								<< negScoreMax_[i*M] << std::endl;
		}

		ofile_mops_logOdds 	<< "positive" << '\t'
							<< "negative" << std::endl;
		for( i = 0; i < posScoreAll_.size(); i++ ){
			ofile_mops_logOdds 	<< std::setprecision(6)
								<< posScoreAll_[i] << '\t'
								<< negScoreAll_[i*M] << std::endl;
		}
	}
}
