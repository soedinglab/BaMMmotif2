#include "FDR.h"

#include <float.h>		// -FLT_MAX

FDR::FDR( Motif* motif ){

	motif_ = motif;

	trainFolds_.resize( Global::cvFold - 1 );
	testFold_.resize( 1 );

	for( int k = 0; k < std::max( Global::modelOrder+2, 4 ); k++ ){	// 4 is for cases when modelOrder < 2
		Y_.push_back( ipow( Alphabet::getSize(), k ) );
	}

	// TODO: only for testing:
	std::cout << "Motif parameter before CV, v[k][y][1]:\n";
	for( int k = 0; k < Global::modelOrder+1; k++ ){
		for( int y = 0; y < Y_[k+1]; y++ ){
			std::cout << motif_->getV()[k][y][1] << '\t';
		}
		std::cout << std::endl;
	}

	int posN = Global::posSequenceSet->getN();
	int negN = posN * Global::mFold;
	int LW1 = Global::posSequenceSet->getMaxL()-motif_->getW()+1;
	int n;

	posScores_ = ( float** )calloc( posN, sizeof( float* ) );
	for( n = 0; n < posN; n++ ){
		posScores_[n] = ( float* )calloc( LW1, sizeof( float ) );
	}
	negScores_ = ( float** )calloc( negN, sizeof( float* ) );
	for( n = 0; n < negN; n++ ){
		negScores_[n] = ( float* )calloc( LW1, sizeof( float ) );
	}
	FP_mops_ = new float[( posN + negN ) * LW1];
	TFP_mops_= new float[( posN + negN ) * LW1];

}

FDR::~FDR(){

	for( int n = 0; n < Global::posSequenceSet->getN(); n++ ){
		free( posScores_[n] );
	}
	free( posScores_ );

	for( int n = 0; n < Global::posSequenceSet->getN() * Global::mFold; n++ ){
		free( negScores_[n] );
	}
	free( negScores_ );

	delete[] FP_mops_;
	delete[] TFP_mops_;
}

void FDR::evaluateMotif(){

	for( int fold = 0; fold < Global::cvFold; fold++ ){

		printf( " ________________________________\n"
				"|                                |\n"
				"|  Cross validation for fold-%d   |\n"
				"|________________________________|\n\n", fold+1 );

		Motif* fullMotif = motif_;

		// TODO: only for testing:
		std::cout << "Motif parameter before each fold of CV, v[k][y][1]:\n";

		for( int k = 0; k < Global::modelOrder+1; k++ ){
			for( int y = 0; y < Y_[k+1]; y++ ){
				std::cout << fullMotif->getV()[k][y][1] << '\t';
			}
			std::cout << std::endl;
		}

		// assign training and test folds
		trainFolds_.clear();
		testFold_.clear();
		for( int f = 0; f < Global::cvFold; f++ ){
			if( f != fold ){
				trainFolds_.push_back( f );
			}
		}
		testFold_.push_back( fold );

		BackgroundModel* trainBgModel = new BackgroundModel( *Global::negSequenceSet,
													Global::bgModelOrder,
													Global::bgModelAlpha,
													Global::interpolateBG,
													Global::posFoldIndices,
													trainFolds_ );

		BackgroundModel* testBgModel = new BackgroundModel( *Global::posSequenceSet,
													Global::bgModelOrder,
													Global::bgModelAlpha,
													Global::interpolateBG,
													Global::posFoldIndices,
													testFold_ );

		EM* em = new EM( fullMotif, trainBgModel, trainFolds_ );
		em->learnMotif();

		long timestamp1 = time( NULL );
		// score positive test sequences
		scoreSequenceSet( fullMotif, trainBgModel, testFold_ );
		fprintf( stdout, "\n--- Runtime for scoring positive test set: %ld seconds ---\n", time( NULL )-timestamp1 );

		long timestamp2 = time( NULL );
		// generate negative sequences
		std::vector<Sequence*> negTestSequenceSet = sampleSequenceSet( testBgModel->getV() );
		fprintf( stdout, "\n--- Runtime for generating negative set: %ld seconds ---\n", time( NULL )-timestamp2 );

		long timestamp3 = time( NULL );
		// score negative sequences
		scoreSequenceSet( fullMotif, trainBgModel, negTestSequenceSet );
		fprintf( stdout, "\n--- Runtime for scoring negative sequence set: %ld seconds ---\n", time( NULL )-timestamp3 );
	}

	long timestamp4 = time( NULL );
	printf( " __________________________________\n"
			"|                                  |\n"
			"|  calculate precision and recall  |\n"
			"|__________________________________|\n\n" );
	calculatePR();
	fprintf( stdout, "\n--- Runtime for calculating precision and recall: %ld seconds ---\n", time( NULL )-timestamp4 );
}

// score Global::posSequenceSet using Global::posFoldIndices and folds
void FDR::scoreSequenceSet( Motif* motif, BackgroundModel* bg, std::vector<int> testFold ){

	int K_motif = Global::modelOrder;
	int K_bg = Global::bgModelOrder;
	int W = motif->getW();
	int i, LW1, k, y, j;
	std::vector<int> testIndices = Global::posFoldIndices[testFold[0]];
	std::vector<Sequence*> posSeqs = Global::posSequenceSet->getSequences();
	float maxS;												// maximal logOddsScore over all positions for each sequence

	for( size_t n = 0; n < testIndices.size(); n++ ){
		LW1 = posSeqs[testIndices[n]]->getL() - W + 1;
		// reset maxS to a minimum negative number
		maxS = -FLT_MAX;
		for( i = 0; i < LW1; i++ ){
			for( j = 0; j < W; j++ ){
				k = std::min( i+j, K_motif );
				y = posSeqs[testIndices[n]]->extractKmer( i+j, k );
				posScores_[testIndices[n]][i] += ( logf( motif->getV()[K_motif][y][j] )
											- logf( bg->getV()[std::min( K_motif, K_bg )][y] ) );
			}

			// take the complete logOddsScores for MOPS model:
			posScoresAll_.push_back( posScores_[testIndices[n]][i] );

			// take the largest logOddsScore for ZOOPS model:
			if( posScores_[testIndices[n]][i] > maxS ){
				maxS = posScores_[testIndices[n]][i];
			}
		}
		posScoresMax_.push_back( maxS );
	}

}

// score sequences from negative sequence set
void FDR::scoreSequenceSet( Motif* motif, BackgroundModel* bg, std::vector<Sequence*> seqSet ){

	int K_motif = Global::modelOrder;
	int K_bg = Global::bgModelOrder;
	int W = motif->getW();
	int i, LW1, k, y, j;
	std::vector<Sequence*> seqs = seqSet;
	float maxS;												// maximal logOddsScore over all positions for each sequence

	for( size_t n = 0; n < seqs.size(); n++ ){
		LW1 = seqs[n]->getL() - W + 1;
		maxS = -FLT_MAX;
		for( i = 0; i < LW1; i++ ){
			for( j = 0; j < W; j++ ){
				k = std::min(i+j, K_motif );
				y = seqs[n]->extractKmer( i+j, k );
				negScores_[n][i] += ( logf( motif->getV()[K_motif][y][j] ) - logf( bg->getV()[std::min( K_motif, K_bg )][y] ) );
			}

			// take the complete logOddsScores for MOPS model:
			negScoresAll_.push_back( negScores_[n][i] );

			// take the largest logOddsScore for ZOOPS model:
			if( negScores_[n][i] > maxS ){
				maxS = negScores_[n][i];
			}
		}
		negScoresMax_.push_back( maxS );
	}
}

//// generate negative sequences based on full positive set
//std::vector<Sequence*> FDR::sampleSequenceSet(){
//
//	int posN = Global::posSequenceSet->getN();
//	int negN = posN * Global::mFold;
//	float** freqs = Global::posSequenceSet->getKmerFrequencies();
//
//	// generate sample sequence set based on the (k+1)-mer frequencies
//	std::vector<Sequence*> sampleSet;
//
//	for( int n = 0; n < negN; n++ ){
//
//		// generate sample sequence
//		int L = Global::posSequenceSet->getMaxL();  		// TODO: better to use the average length of all the sequences
//		uint8_t* sequence = ( uint8_t* )calloc( L, sizeof( uint8_t ) );
//		std::string header = "> sample sequence";
//
//		int i, k, yk;
//		double f, random;
//		uint8_t a;
//		int K = std::min( 5, Global::modelOrder );			// Fix the kmer frequencies to order 2
//
//		// get a random number for the first nucleotide
//		random = ( double )rand() / ( double )RAND_MAX;
//		f = 0.0f;
//		for( a = 1; a <= Y_[1]; a++ ){
//			f += freqs[0][a-1];
//			if( random <= f ){
//				sequence[0] = a;
//				break;
//			}
//			if( sequence[0] == 0 )	sequence[0] = a;		// Trick: this is to solve the numerical problem
//		}
//
//		for( i = 1; i < L; i++ ){
//			random = ( double )rand() / ( double )RAND_MAX;	// get another random double number
//			// calculate y of K-mer
//			yk = 0;
//			for( k = std::min( i, K ); k > 0; k-- ){
//				yk += ( sequence[i-k] - 1 ) * Y_[k];
//			}
//
//			// assign a nucleotide based on K-mer frequency
//			f = 0.0f;
//			for( a = 1; a <= Y_[1]; a++ ){
//				f += freqs[std::min( i, K )][yk+a-1];
//				if( random <= f ){
//					sequence[i] = a;
//					break;
//				}
//				if( sequence[i] == 0 )	sequence[i] = a;	// Trick: this is to solve the numerical problem
//			}
//		}
//
//		Sequence* sampleSequence = new Sequence( sequence, L, header, Y_, Global::revcomp );
//		sampleSet.push_back( sampleSequence );
//	}
//
//	return sampleSet;
//
//}


// generate negative sequences based on each test set
std::vector<Sequence*> FDR::sampleSequenceSet( float** v ){

	int posN = Global::posSequenceSet->getN();
	int negN = posN * Global::mFold / Global::cvFold;	// TODO: default Global::mFold / Global::cvFold to an integer

	// generate sample sequence set based on the (k+1)-mer frequencies
	std::vector<Sequence*> sampleSet;
	for( int i = 0; i < negN; i++ ){
		sampleSet.push_back( sampleSequence( v ) );
	}
	return sampleSet;
}

// generate sample sequence based on k-mer probabilities
Sequence* FDR::sampleSequence( float** v ){

	int L = Global::posSequenceSet->getMaxL();  		// TODO: better to use the average length of all the sequences
	uint8_t* sequence = ( uint8_t* )calloc( L, sizeof( uint8_t ) );
	std::string header = "sample sequence";

	int i, k, yk;
	double f, random;
	uint8_t a;
	int K = Global::modelOrder;

	// get a random number for the first nucleotide
	random = ( double )rand() / ( double )RAND_MAX;
	f = 0.0f;
	for( a = 1; a <= Y_[1]; a++ ){
		f += v[0][a-1];
		if( random <= f ){
			sequence[0] = a;
			break;
		}
		if( sequence[0] == 0 )	sequence[0] = a;		// Trick: this is to solve the numerical problem
	}

	for( i = 1; i < L; i++ ){
		random = ( double )rand() / ( double )RAND_MAX;	// get another random double number
		// calculate y of K-mer
		yk = 0;
		for( k = std::min( i, K ); k > 0; k-- ){
			yk += ( sequence[i-k] - 1 ) * Y_[k];
		}

		// assign a nucleotide based on K-mer frequency
		f = 0.0f;
		for( a = 1; a <= Y_[1]; a++ ){
			f += v[std::min( i, K )][yk+a-1];
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
	// for ZOOPS model:
	std::sort( posScoresMax_.begin(), posScoresMax_.end(), std::greater<float>() );
	std::sort( negScoresMax_.begin(), negScoresMax_.end(), std::greater<float>() );
	// for MOPS model:
	std::sort( posScoresAll_.begin(), posScoresAll_.end(), std::greater<float>() );
	std::sort( negScoresAll_.begin(), negScoresAll_.end(), std::greater<float>() );

	// Rank and score these log odds score values
	size_t size_max = posScoresMax_.size() + negScoresMax_.size();
	size_t size_all = posScoresAll_.size() + negScoresAll_.size();
	size_t i;
	size_t i_posM = 0, i_negM = 0;						// index for arrays storing the max log odds scores
	size_t i_posA = 0, i_negA = 0;						// index for arrays storing the complete log odds scores
	size_t N = Global::posSequenceSet->getN();

	// For "many occurrence per sequence" (MOPS) model
	float max_diff = 0.0f;								// maximal difference between TFP and FP
	size_t idx_max = 0;									// index when the difference between TFP and FP reaches maximum

	for( i = 0; i < size_all; i++ ){
		if( posScoresAll_[i_posA] > negScoresAll_[i_negA] ||
				i_posA == 0 || i_negA == negScoresAll_.size() ){
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
			idx_max = 20 * N;
			break;
		}
	}

	for( i = 0; i < idx_max; i++ ){
		P_mops_.push_back( 1 - FP_mops_[i] / TFP_mops_[i] );
		R_mops_.push_back( ( TFP_mops_[i] - FP_mops_[i] ) / max_diff );
	}

	// For "zero or one occurrence per sequence" (ZOOPS) model
	for( i = 0; i < size_max; i++ ){
		if( posScoresMax_[i_posM] > negScoresMax_[i_negM] ||
				i_posM == 0 || i_negM == negScoresMax_.size() ){
			i_posM++;
		} else {
			i_negM++;
		}

		P_zoops_.push_back( static_cast<float>( i_posM ) / static_cast<float>( i+1 ) );  // i_posM is the true positive value
		R_zoops_.push_back( static_cast<float>( i_posM ) / static_cast<float>( N ) );
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
	 */

	std::string opath = std::string( Global::outputDirectory ) + '/'
			+ std::string( Global::posSequenceBasename );
	std::string opath_zoops_p = opath + ".zoops.precision";
	std::string opath_zoops_r = opath + ".zoops.recall";
	std::string opath_mops_p = opath + ".mops.precision";
	std::string opath_mops_r = opath + ".mops.recall";
	std::string opath_mops_fp = opath + ".mops.fp";
	std::string opath_mops_tfp = opath + ".mops.tfp";

	std::ofstream ofile_zoops_p( opath_zoops_p );
	std::ofstream ofile_zoops_r( opath_zoops_r );
	std::ofstream ofile_mops_p( opath_mops_p );
	std::ofstream ofile_mops_r( opath_mops_r );
	std::ofstream ofile_mops_fp( opath_mops_fp );
	std::ofstream ofile_mops_tfp( opath_mops_tfp );

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
}
