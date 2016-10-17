#include "FDR.h"

#include <float.h>		// -FLT_MAX

FDR::FDR( Motif* motif ){

	motif_ = motif;

	trainsetFolds_.resize( Global::cvFold - 1 );

	for( int k = 0; k < std::max( Global::modelOrder+2 , 4 ); k++ ){	// 4 is for cases when modelOrder < 2
		Y_.push_back( ipow( Alphabet::getSize(), k ) );
	}

	int posN = Global::posSequenceSet->getN();
	int LW1 = Global::posSequenceSet->getMaxL()-motif_->getW()+1;		// TODO: using maxL can be a waste of memory

	FP_mops_ = new float[posN * ( Global::mFold+1 ) * LW1];
	TFP_mops_= new float[posN * ( Global::mFold+1 ) * LW1];

	testsetV_ = ( float** )calloc( 3, sizeof( float* ) );			// fix to trimer frequencies
	testsetN_ = ( int** )calloc( 3, sizeof( int* ) );
	for( int k = 0; k < 3; k++ ){
		testsetV_[k] = ( float* )calloc( Y_[k+1], sizeof( float ) );
		testsetN_[k] = ( int* )calloc( Y_[k+1], sizeof( int ) );
	}

}

FDR::~FDR(){

	delete[] FP_mops_;
	delete[] TFP_mops_;

	for( int k = 0; k < Global::modelOrder+1; k++ ){
		free( testsetV_[k] );
		free( testsetN_[k] );
	}
	free( testsetV_ );
	free( testsetN_ );
}

void FDR::evaluateMotif(){

	std::vector<Sequence*> posSeqs = Global::posSequenceSet->getSequences();

	for( int fold = 0; fold < Global::cvFold; fold++ ){

		printf( " ________________________________\n"
				"|                                |\n"
				"|  Cross validation for fold-%d   |\n"
				"|________________________________|\n\n", fold+1 );

		Motif* motif = new Motif( *motif_ );			// deep copy motif

		// assign training folds
		trainsetFolds_.clear();
		for( int f = 0; f < Global::cvFold; f++ ){
			if( f != fold ){
				trainsetFolds_.push_back( f );
			}
		}

		// draw sequences for each test set
		std::vector<Sequence*> testSet;
		for( size_t i = 0; i < Global::posFoldIndices[fold].size(); i++ ){
			testSet.push_back( posSeqs[Global::posFoldIndices[fold][i]]);
		}

		// obtain background model for each training set
		BackgroundModel* bgModel = new BackgroundModel( *Global::negSequenceSet,
													Global::bgModelOrder,
													Global::bgModelAlpha,
													Global::interpolateBG,
													Global::posFoldIndices,
													trainsetFolds_ );

		EM* em = new EM( motif, bgModel, trainsetFolds_ );
		em->learnMotif();

		clock_t t0 = clock();
		// score positive test sequences
		posSetScores_= scoreSequenceSet( motif, bgModel, testSet );
		fprintf( stdout, "\n--- Runtime for scoring positive test set: %.4f seconds ---\n", ( ( float )( clock() - t0 ) ) / CLOCKS_PER_SEC );

		if( !Global::bgSequenceFilename ){
			t0 = clock();
			// generate negative sequences
			std::vector<Sequence*> negTestSequenceSet = sampleSequenceSet( testSet );
			fprintf( stdout, "\n--- Runtime for generating negative set: %.4f seconds ---\n", ( ( float )( clock() - t0 ) ) / CLOCKS_PER_SEC );

			t0 = clock();
			// score negative sequences
			negSetScores_ = scoreSequenceSet( motif, bgModel, negTestSequenceSet );
			fprintf( stdout, "\n--- Runtime for scoring negative sequence set: %.4f seconds ---\n", ( ( float )( clock() - t0 ) ) / CLOCKS_PER_SEC );
		}
	}

	clock_t t0 = clock();
	if( Global::bgSequenceFilename ){
		Motif* motif = new Motif( *motif_ );
		BackgroundModel* bgModel = new BackgroundModel( *Global::negSequenceSet,
													Global::bgModelOrder,
													Global::bgModelAlpha );
		negSetScores_ = scoreSequenceSet( motif, bgModel, Global::bgSequenceSet->getSequences() );
		fprintf( stdout, "\n--- Runtime for scoring negative sequence set: %.4f seconds ---\n", ( ( float )( clock() - t0 ) ) / CLOCKS_PER_SEC );
	}

	t0 = clock();
	printf( " __________________________________\n"
			"|                                  |\n"
			"|  calculate precision and recall  |\n"
			"|__________________________________|\n\n" );
	calculatePR();
	fprintf( stdout, "\n--- Runtime for calculating precision and recall: %.4f seconds ---\n", ( ( float )( clock() - t0 ) ) / CLOCKS_PER_SEC );
}

// score sequences for both positive and negative sequence sets
std::vector<std::vector<float>> FDR::scoreSequenceSet( Motif* motif, BackgroundModel* bg, std::vector<Sequence*> seqSet ){

	std::vector<std::vector<float>> scores( 2 );			// scores[0]: store the log odds scores at all positions of each sequence
															// scores[1]: store maximal log odds scores of each sequence
	int K_motif = Global::modelOrder;
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
				int y = seqSet[n]->extractKmer( i+j, std::min(i+j, K_motif ) );
				int y_bg = y % Y_[K_bg+1];
				logOdds[i] += ( logf( motif->getV()[K_motif][y][j] ) - logf( bg->getV()[std::min( K_motif, K_bg )][y_bg] ) );
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
std::vector<Sequence*> FDR::sampleSequenceSet( std::vector<Sequence*> seqs ){

	std::vector<Sequence*> sampleSet;

	calculateTrimerV( seqs );

	for( size_t i = 0; i < seqs.size(); i++ ){
		int L = seqs[i]->getL();
		for( int n = 0; n < Global::mFold; n++ ){
			sampleSet.push_back( sampleSequence( L, testsetV_ ) );
		}
	}
	return sampleSet;
}

// generate sample sequence based on trimer conditional probabilities
Sequence* FDR::sampleSequence( int L, float** v ){

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
		for( int k = std::min( i, 2 ); k > 0; k-- ){
			yk += ( sequence[i-k] - 1 ) * Y_[k];
		}

		// assign a nucleotide based on K-mer frequency
		f = 0.0f;
		for( uint8_t a = 1; a <= Y_[1]; a++ ){
			f += v[std::min( i, 2 )][yk+a-1];
			if( random <= f ){
				sequence[i] = a;
				break;
			}
			if( sequence[i] == 0 )	sequence[i] = a;	// Trick: this is to solve the numerical problem
		}
	}

	Sequence* sampleSequence = new Sequence( sequence, L, header, Y_, Global::revcomp );

	return sampleSequence;
}

void FDR::calculatePR(){

	// Sort log odds scores from large to small
	// for MOPS model:
	std::sort( posSetScores_[0].begin(), posSetScores_[0].end(), std::greater<float>() );
	std::sort( negSetScores_[0].begin(), negSetScores_[0].end(), std::greater<float>() );
	// for ZOOPS model:
	std::sort( posSetScores_[1].begin(), posSetScores_[1].end(), std::greater<float>() );
	std::sort( negSetScores_[1].begin(), negSetScores_[1].end(), std::greater<float>() );

	// Rank and score these log odds score values
	size_t size_all = posSetScores_[0].size() + negSetScores_[0].size();
	size_t size_max = posSetScores_[1].size() + negSetScores_[1].size();
	size_t i_posM = 0, i_negM = 0;						// index for arrays storing the max log odds scores
	size_t i_posA = 0, i_negA = 0;						// index for arrays storing the complete log odds scores
	size_t N = Global::posSequenceSet->getN() / Global::cvFold;

	// For "many occurrence per sequence" (MOPS) model
	float max_diff = 0.0f;								// maximal difference between TFP and FP
	size_t idx_max = 0;									// index when the difference between TFP and FP reaches maximum

	int posN = Global::posSequenceSet->getN();
	int LW1 = Global::posSequenceSet->getMaxL()-motif_->getW()+1;
	FP_mops_ = new float[posN * ( Global::mFold+1 ) * LW1];
	TFP_mops_= new float[posN * ( Global::mFold+1 ) * LW1];

	for( size_t i = 0; i < size_all; i++ ){
		if( posSetScores_[0][i_posA] > negSetScores_[0][i_negA] ||
				i_posA == 0 || i_negA == negSetScores_[0].size() ){
			i_posA++;
		} else {
			i_negA++;
		}

		TFP_mops_[i] = static_cast<float>( i_posA );
		FP_mops_[i] = static_cast<float>( i_negA ) / static_cast<float>( Global::mFold );

		if( max_diff == TFP_mops_[i] - FP_mops_[i] ) {	// stop when recall reaches 1
			idx_max = i;
		}
		if( max_diff < TFP_mops_[i] - FP_mops_[i] ){
			max_diff = TFP_mops_[i] - FP_mops_[i];
		}
		if( idx_max > 20 * N ){
			break;
		}
	}

	for( size_t i = 0; i < idx_max; i++ ){
		P_mops_.push_back( 1 - FP_mops_[i] / TFP_mops_[i] );
		R_mops_.push_back( ( TFP_mops_[i] - FP_mops_[i] ) / max_diff );
	}

	// For "zero or one occurrence per sequence" (ZOOPS) model
	for( size_t i = 0; i < size_max; i++ ){
		if( posSetScores_[1][i_posM] > negSetScores_[1][i_negM] ||
				i_posM == 0 || i_negM == negSetScores_[1].size() ){
			i_posM++;
		} else {
			i_negM++;
		}

		P_zoops_.push_back( static_cast<float>( i_posM ) / static_cast<float>( i+1 ) );  // i_posM = TP
		R_zoops_.push_back( static_cast<float>( i_posM ) / static_cast<float>( N ) );
	}

}

void FDR::calculateTrimerV( std::vector<Sequence*> seqs ){

	// reset counts for trimers
	for( int k = 0; k < 3; k++ ){
		for( int y = 0; y < Y_[k+1]; y++ ){
			testsetN_[k][y] = 0;
		}
	}

	// count trimers
	for( size_t i = 0; i < seqs.size(); i++ ){
		int L = seqs[i]->getL();
		for( int k = 0; k < 3; k++ ){
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

	// calculate conditional probabilities for kmers till 3rd order
	int normFactor = 0;
	for( int y = 0; y < Y_[1]; y++ )	normFactor += testsetN_[0][y];
	for( int y = 0; y < Y_[1]; y++ )	testsetV_[0][y] = ( float )testsetN_[0][y] / ( float )normFactor;
	for( int k = 1; k < 3; k++ ){
		for( int y = 0; y < Y_[k+1]; y++ ){
			int yk = y / Y_[1];
			testsetV_[k][y] = ( float )testsetN_[k][y] / ( float )testsetN_[k-1][yk];
		}
	}
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
	 * (7) posSequenceBasename.zoops.posLogOdds
	 * (8) posSequenceBasename.mops.posLogOdds
	 * (9) posSequenceBasename.zoops.negLogOdds
	 * (10)posSequenceBasename.mops.negLogOdds
	 */

	std::string opath = std::string( Global::outputDirectory ) + '/'
			+ std::string( Global::posSequenceBasename );
	std::string opath_zoops_p = opath + ".zoops.precision";
	std::string opath_zoops_r = opath + ".zoops.recall";
	std::string opath_mops_p = opath + ".mops.precision";
	std::string opath_mops_r = opath + ".mops.recall";
	std::string opath_mops_fp = opath + ".mops.fp";
	std::string opath_mops_tfp = opath + ".mops.tfp";
	std::string opath_zoops_posScore = opath + ".zoops.posLogOdds";
	std::string opath_mops_posScore = opath + ".mops.posLogOdds";
	std::string opath_zoops_negScore = opath + ".zoops.negLogOdds";
	std::string opath_mops_negScore = opath + ".mops.negLogOdds";

	std::ofstream ofile_zoops_p( opath_zoops_p );
	std::ofstream ofile_zoops_r( opath_zoops_r );
	std::ofstream ofile_mops_p( opath_mops_p );
	std::ofstream ofile_mops_r( opath_mops_r );
	std::ofstream ofile_mops_fp( opath_mops_fp );
	std::ofstream ofile_mops_tfp( opath_mops_tfp );
	std::ofstream ofile_zoops_posScore( opath_zoops_posScore );
	std::ofstream ofile_mops_posScore( opath_mops_posScore );
	std::ofstream ofile_zoops_negScore( opath_zoops_negScore );
	std::ofstream ofile_mops_negScore( opath_mops_negScore );

	size_t i;
	size_t N = Global::posSequenceSet->getN();

	for( i = 0; i < P_zoops_.size(); i++ ){
		ofile_zoops_p << P_zoops_[i] << std::endl;
		ofile_zoops_r << R_zoops_[i] << std::endl;
	}

	for( i = 0; i < P_mops_.size(); i++ ){
		ofile_mops_p << P_mops_[i] << std::endl;
		ofile_mops_r << R_mops_[i] << std::endl;
	}

	for( i = 0; i < 20 * N; i++ ){
		ofile_mops_fp << FP_mops_[i] << std::endl;
		ofile_mops_tfp << TFP_mops_[i] << std::endl;
	}

	for( i = 0; i < posSetScores_[1].size(); i++ ){
		ofile_zoops_posScore << posSetScores_[1][i] << std::endl;
	}

	for( i = 0; i < posSetScores_[0].size(); i++ ){
		ofile_mops_posScore << posSetScores_[0][i] << std::endl;
	}

	for( i = 0; i < negSetScores_[1].size(); i++ ){
		ofile_zoops_negScore << negSetScores_[1][i] << std::endl;
	}

	for( i = 0; i < negSetScores_[0].size(); i++ ){
		ofile_mops_negScore << negSetScores_[0][i] << std::endl;
	}
}
