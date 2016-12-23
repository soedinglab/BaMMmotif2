#include "FDR.h"

#include <float.h>		// -FLT_MAX

FDR::FDR( Motif* motif ){

	motif_ = motif;

	for( int k = 0; k < std::max( Global::modelOrder+2 , Global::bgModelOrder+2 ); k++ ){
		Y_.push_back( ipow( Alphabet::getSize(), k ) );
	}

	testsetV_ = ( float** )calloc( Global::sOrder+1, sizeof( float* ) );
	testsetN_ = ( int** )calloc( Global::sOrder+1, sizeof( int* ) );
	for( int k = 0; k < Global::sOrder+1; k++ ){
		testsetV_[k] = ( float* )calloc( Y_[k+1], sizeof( float ) );
		testsetN_[k] = ( int* )calloc( Y_[k+1], sizeof( int ) );
	}
}

FDR::~FDR(){

	for( int k = 0; k < Global::sOrder+1; k++ ){
		free( testsetV_[k] );
		free( testsetN_[k] );
	}
	free( testsetV_ );
	free( testsetN_ );
}

void FDR::evaluateMotif(){

	std::vector<Sequence*> posSeqs = Global::posSequenceSet->getSequences();

	std::vector<std::vector<float>> scores;
	for( int fold = 0; fold < Global::cvFold; fold++ ){

		printf( " ________________________________\n"
				"|                                |\n"
				"|  Cross validation for fold-%d   |\n"
				"|________________________________|\n\n", fold+1 );

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
		ModelLearning model( motif, bgModel, trainsetFolds );
		if( Global::EM ){
			model.EMlearning();
		} else if( Global::CGS ){
			model.GibbsSampling();
		}

		// score positive test sequences on optimized motif
		scores = scoreSequenceSet( motif, bgModel, testSet );
		posScoreAll_.insert( std::end( posScoreAll_ ), std::begin( scores[0] ), std::end( scores[0] ) );
		posScoreMax_.insert( std::end( posScoreMax_ ), std::begin( scores[1] ), std::end( scores[1] ) );

		std::vector<std::unique_ptr<Sequence>> negSet;

		// generate negative sequence set
		negSet = sampleSequenceSet( testSet );

		// score negative sequence set
		scores = scoreSequenceSet( motif, bgModel, negSet );
		negScoreAll_.insert( std::end( negScoreAll_ ), std::begin( scores[0] ), std::end( scores[0] ) );
		negScoreMax_.insert( std::end( negScoreMax_ ), std::begin( scores[1] ), std::end( scores[1] ) );

		delete motif;
		delete bgModel;
	}

	printf( " __________________________________\n"
			"|                                  |\n"
			"|  calculate precision and recall  |\n"
			"|__________________________________|\n\n" );
	calculatePR();
}

// score sequences for both positive and negative sequence sets
std::vector<std::vector<float>> FDR::scoreSequenceSet( Motif* motif, BackgroundModel* bg, const std::vector<std::unique_ptr<Sequence>> & seqSet ){

	std::vector<std::vector<float>> scores( 2 );			// scores[0]: store the log odds scores at all positions of each sequence
															// scores[1]: store maximal log odds scores of each sequence
	int K = Global::modelOrder;
	int K_bg = Global::bgModelOrder;
	int W = motif->getW();
	float maxScore;											// maximal logOddsScore over all positions for each sequence

	for( size_t n = 0; n < seqSet.size(); n++ ){
		int LW1 = seqSet[n]->getL() - W + 1;
		maxScore = -FLT_MAX;
		std::vector<float> logOdds( LW1 );					// calculate log odds scores for all the positions
		for( int i = 0; i < LW1; i++ ){
			logOdds[i] = 0.0f;
			for( int j = 0; j < W; j++ ){
				int y = seqSet[n]->extractKmer( i+j, std::min(i+j, K ) );
				int y_bg = y % Y_[K_bg+1];
				if( y >= 0 ){
					logOdds[i] += ( logf( motif->getV()[K][y][j] ) - logf( bg->getV()[std::min( K, K_bg )][y_bg] ) );
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

// score sequences for both positive and negative sequence sets
std::vector<std::vector<float>> FDR::scoreSequenceSet( Motif* motif, BackgroundModel* bg, std::vector<Sequence*> seqSet ){

	std::vector<std::vector<float>> scores( 2 );			// scores[0]: store the log odds scores at all positions of each sequence
															// scores[1]: store maximal log odds scores of each sequence
	int K = Global::modelOrder;
	int K_bg = Global::bgModelOrder;
	int W = motif->getW();
	float maxScore;											// maximal logOddsScore over all positions for each sequence

	for( size_t n = 0; n < seqSet.size(); n++ ){
		int LW1 = seqSet[n]->getL() - W + 1;
		maxScore = -FLT_MAX;
		std::vector<float> logOdds( LW1 );					// calculate log odds scores for all the positions
		for( int i = 0; i < LW1; i++ ){
			logOdds[i] = 0.0f;
			for( int j = 0; j < W; j++ ){
				int y = seqSet[n]->extractKmer( i+j, std::min(i+j, K ) );
				int y_bg = y % Y_[K_bg+1];
				if( y >= 0 ){
					logOdds[i] += ( logf( motif->getV()[K][y][j] ) - logf( bg->getV()[std::min( K, K_bg )][y_bg] ) );
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
// generate negative sequences based on each test set
std::vector<std::unique_ptr<Sequence>> FDR::sampleSequenceSet( std::vector<Sequence*> seqs ){

	std::vector<std::unique_ptr<Sequence>> sampleSet;

	calcKmerFreq( seqs );

	for( size_t i = 0; i < seqs.size(); i++ ){
		int L = seqs[i]->getL();
		for( int n = 0; n < Global::mFold; n++ ){
			sampleSet.push_back( sampleSequence( L, testsetV_ ) );
		}
	}
	return sampleSet;
}

// generate sample sequence based on trimer conditional probabilities
std::unique_ptr<Sequence> FDR::sampleSequence( int L, float** v ){

	uint8_t* sequence = ( uint8_t* )calloc( L, sizeof( uint8_t ) );
	std::string header = "sample sequence";

	// get a random number for the first nucleotide
	double random = ( double )rand() / ( double )RAND_MAX;
	double f = 0.0f;
	for( uint8_t a = 1; a <= Y_[1]; a++ ){
		f += v[0][a-1];
		if( random <= f ){
			sequence[0] = a;
			break;
		}
		if( sequence[0] == 0 )	sequence[0] = a;		// Trick: this is to solve the numerical problem
	}

	for( int i = 1; i < L; i++ ){
		random = ( double )rand() / ( double )RAND_MAX;	// get another random double number
		// calculate y of K-mer
		int yk = 0;
		for( int k = std::min( i, Global::sOrder ); k > 0; k-- ){
			yk += ( sequence[i-k] - 1 ) * Y_[k];
		}

		// assign a nucleotide based on K-mer frequency
		f = 0.0f;
		for( uint8_t a = 1; a <= Y_[1]; a++ ){
			f += v[std::min( i, Global::sOrder )][yk+a-1];
			if( random <= f ){
				sequence[i] = a;
				break;
			}
			if( sequence[i] == 0 )	sequence[i] = a;	// Trick: this is to solve the numerical problem
		}
	}
	std::unique_ptr<Sequence> seq( new Sequence( sequence, L, header, Y_, Global::revcomp ) );
	free( sequence );
	return seq;

}

void FDR::calculatePR(){

	int posN = Global::posSequenceSet->getN();
	int negN = posN * Global::mFold;

	// for MOPS model:
	// Sort log odds scores from large to small
	std::sort( posScoreAll_.begin(), posScoreAll_.end(), std::greater<float>() );
	std::sort( negScoreAll_.begin(), negScoreAll_.end(), std::greater<float>() );
	// Rank and score these log odds score values
	int idx_posAll = 0;
	int idx_negAll = 0;									// index for arrays storing the complete log odds scores
	float max_diff = 0.0f;								// maximal difference between TFP and FP
	int idx_max = posN + negN;							// index when the difference between TFP and FP reaches maximum
														// set the initial cutoff as posN+negN

	for( int i = 0; i < posN + negN; i++ ){
		if( posScoreAll_[idx_posAll] > negScoreAll_[idx_negAll] ||
				idx_posAll == 0 || idx_negAll == posN+negN-1 ){
			idx_posAll++;
		} else {
			idx_negAll++;
		}

		TFP_MOPS_.push_back( ( float )idx_posAll );
		FP_MOPS_.push_back( ( float )idx_negAll / ( float )Global::mFold );

		if( max_diff == TFP_MOPS_[i] - FP_MOPS_[i] ){	// stop when recall reaches 1
			idx_max = i;
			break;
		}
		if( max_diff < TFP_MOPS_[i] - FP_MOPS_[i] ){
			max_diff = TFP_MOPS_[i] - FP_MOPS_[i];
		}
	}

	for( int i = 0; i < idx_max; i++ ){
		Pre_MOPS_.push_back( 1.0f - FP_MOPS_[i] / TFP_MOPS_[i] );
		Rec_MOPS_.push_back( ( TFP_MOPS_[i] - FP_MOPS_[i] ) / max_diff );
		if( ( TFP_MOPS_[i] - FP_MOPS_[i] ) / max_diff <= 0.5f ){
			pre_half_ = Pre_MOPS_[i];
		}
	}

	// for ZOOPS model:
	// Sort log odds scores from large to small
	std::sort( posScoreMax_.begin(), posScoreMax_.end(), std::greater<float>() );
	std::sort( negScoreMax_.begin(), negScoreMax_.end(), std::greater<float>() );

	// Rank and score these log odds score values
	int idx_posMax = 0;
	int idx_negMax = 0;									// index for arrays storing the max log odds scores

	for( int i = 0; i < posN + negN; i++ ){
		if( posScoreMax_[idx_posMax] > negScoreMax_[idx_negMax] ||
				idx_posMax == 0 || idx_negMax == posN+negN-1 ){
			idx_posMax++;
		} else {
			idx_negMax++;
		}

		Pre_ZOOPS_.push_back( ( float )idx_posMax / ( float )( i+1 ) );  // i_posM = TP
		Rec_ZOOPS_.push_back( ( float )idx_posMax / ( float )posN );
	}

}

void FDR::calcKmerFreq( std::vector<Sequence*> seqs ){

	// reset counts for k-mers
	for( int k = 0; k < Global::sOrder+1; k++ ){
		for( int y = 0; y < Y_[k+1]; y++ ){
			testsetN_[k][y] = 0;
		}
	}

	// count k-mers
	for( size_t i = 0; i < seqs.size(); i++ ){
		int L = seqs[i]->getL();
		for( int k = 0; k < Global::sOrder+1; k++ ){
			// loop over sequence positions
			for( int j = k; j < L; j++ ){
				// extract (k+1)-mer
				int y = seqs[i]->extractKmer( j, k );
				// skip non-defined alphabet letters
				if( y >= 0 ){
					// count (k+1)mer
					testsetN_[k][y]++;
				}
			}
		}
	}

	// calculate frequencies
	int normFactor = 0;
	for( int y = 0; y < Y_[1]; y++ )	normFactor += testsetN_[0][y];
	for( int y = 0; y < Y_[1]; y++ )	testsetV_[0][y] = ( float )testsetN_[0][y] / ( float )normFactor;
	for( int k = 1; k < Global::sOrder+1; k++ ){
		for( int y = 0; y < Y_[k+1]; y++ ){
			int yk = y / Y_[1];
			testsetV_[k][y] = ( float )testsetN_[k][y] / ( float )testsetN_[k-1][yk];
		}
	}
}

float FDR::getPre_half(){
	return pre_half_;
}

void FDR::print(){

}

void FDR::write(){

	/**
	 * save FDR parameters in six flat files:
	 * (1) posSequenceBasename.zoops.precision:	precision values for ZOOPS model
	 * (2) posSequenceBasename.zoops.recall:	recall values for ZOOPS model
	 * (3) posSequenceBasename.mops.precision:	precision values for MOPS model
	 * (4) posSequenceBasename.mops.recall:		recall values for MOPS model
	 * (5) posSequenceBasename.mops.fp:			false positives for MOPS model
	 * (6) posSequenceBasename.mops.tfp:		true and false positive values for MOPS model
	 */

	std::string opath = std::string( Global::outputDirectory ) + '/'
			+ Global::posSequenceBasename;
	std::string opath_zoops_p = opath + ".zoops.precision";
	std::string opath_zoops_r = opath + ".zoops.recall";
	std::string opath_mops_p = opath + ".mops.precision";
	std::string opath_mops_r = opath + ".mops.recall";
//	std::string opath_mops_fp = opath + ".mops.fp";
//	std::string opath_mops_tfp = opath + ".mops.tfp";

	std::ofstream ofile_zoops_p( opath_zoops_p );
	std::ofstream ofile_zoops_r( opath_zoops_r );
	std::ofstream ofile_mops_p( opath_mops_p );
	std::ofstream ofile_mops_r( opath_mops_r );
//	std::ofstream ofile_mops_fp( opath_mops_fp );
//	std::ofstream ofile_mops_tfp( opath_mops_tfp );

	size_t i;

	for( i = 0; i < Pre_ZOOPS_.size(); i++ ){
		ofile_zoops_p << Pre_ZOOPS_[i] << std::endl;
		ofile_zoops_r << Rec_ZOOPS_[i] << std::endl;
	}

	for( i = 0; i < Pre_MOPS_.size(); i++ ){
		ofile_mops_p << Pre_MOPS_[i] << std::endl;
		ofile_mops_r << Rec_MOPS_[i] << std::endl;
//		ofile_mops_fp << FP_MOPS_[i] << std::endl;
//		ofile_mops_tfp << TFP_MOPS_[i] << std::endl;
	}


}

void FDR::writeLogOdds(){

	/**
	 * save log odds scores into four files before sorting them
	 * (1) posSequenceBasename.zoops.posLogOdds
	 * (2) posSequenceBasename.mops.posLogOdds
	 * (3) posSequenceBasename.zoops.negLogOdds
	 * (4) posSequenceBasename.mops.negLogOdds
	 */

	std::string opath = std::string( Global::outputDirectory ) + '/'
			+ Global::posSequenceBasename;

	std::string opath_zoops_posScore = opath + ".zoops.posLogOdds";
	std::string opath_mops_posScore = opath + ".mops.posLogOdds";
	std::string opath_zoops_negScore = opath + ".zoops.negLogOdds";
	std::string opath_mops_negScore = opath + ".mops.negLogOdds";
	std::ofstream ofile_zoops_posScore( opath_zoops_posScore );
	std::ofstream ofile_mops_posScore( opath_mops_posScore );
	std::ofstream ofile_zoops_negScore( opath_zoops_negScore );
	std::ofstream ofile_mops_negScore( opath_mops_negScore );

	size_t i;

	for( i = 0; i < posScoreMax_.size(); i++ ){
		ofile_zoops_posScore << posScoreMax_[i] << std::endl;
	}

	for( i = 0; i < posScoreAll_.size(); i+=Global::mFold ){
		ofile_mops_posScore << posScoreAll_[i] << std::endl;
	}

	for( i = 0; i < negScoreMax_.size(); i+=Global::mFold ){
		ofile_zoops_negScore << negScoreMax_[i] << std::endl;
	}

	for( i = 0; i < negScoreAll_.size(); i+=Global::mFold*Global::mFold ){
		ofile_mops_negScore << negScoreAll_[i] << std::endl;
	}
}
