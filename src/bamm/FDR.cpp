#include "FDR.h"

#include <float.h>		// -FLT_MAX

FDR::FDR( Motif* motif ){

	motif_ = motif;

	trainFolds_.resize( Global::cvFold - 1 );
	testFold_ = 0;

	int posN = Global::posSequenceSet->getN();
	int negN = posN * Global::mFold;
	int LW1 = Global::posSequenceSet->getMaxL()-motif_->getW()+1;
	int k, n;

	freqs_ = ( float** )calloc( Global::modelOrder, sizeof( float* ) );
	for( k = 0; k < Global::modelOrder + 1; k++)
		freqs_[k] = ( float* )calloc( ipow( Alphabet::getSize(), k+1 ), sizeof( float ) );
	posS_ = ( float** )calloc( posN, sizeof( float* ) );
	for( n = 0; n < posN; n++ )
		posS_[n] = ( float* )calloc( LW1, sizeof( float ) );

	negS_ = ( float** )calloc( negN, sizeof( float* ) );
	for( n = 0; n < negN; n++ )
		negS_[n] = ( float* )calloc( LW1, sizeof( float ) );

	FP_mops_ = new float[( posN + negN ) * LW1];
	TFP_mops_= new float[( posN + negN ) * LW1];

}

FDR::~FDR(){

	for( int k = 0;  k < Global::modelOrder; k++ ){
		free( freqs_[k] );
	}
	free( freqs_ );
	for( int n = 0; n < Global::posSequenceSet->getN(); n++ ){
		free( posS_[n] );
	}
	free( posS_ );
	for( int n = 0; n < Global::posSequenceSet->getN() * Global::mFold; n++ ){
		free( negS_[n] );
	}
	free( negS_ );

	delete[] FP_mops_;
	delete[] TFP_mops_;
}

void FDR::evaluateMotif(){

	for( int fold = 0; fold < Global::cvFold; fold++ ){

		// assign training and test folds
		trainFolds_.clear();

		for( int f=0; f < Global::cvFold; f++ ){
			if( f != fold ){
				trainFolds_.push_back( f );
			}
		}
		testFold_ =  fold;

		BackgroundModel* bg = new BackgroundModel( *Global::negSequenceSet,
													Global::bgModelOrder,
													Global::bgModelAlpha,
													Global::interpolateBG,
													Global::posFoldIndices,
													trainFolds_ );

		EM* em = new EM( motif_, bg, trainFolds_ );
		em->learnMotif();
		// TODO: to be deleted, only for testing
//		em->write();

		// score positive test sequences
		printf( " _____________________________\n"
				"|                             |\n"
				"|  score test set for fold-%d  |\n"
				"|_____________________________|\n\n", fold );
		scoreSequenceSet( motif_, bg, testFold_ );

		// generate negative sequences
		printf( " __________________________________\n"
				"|                                  |\n"
				"|  create negative set for fold-%d  |\n"
				"|__________________________________|\n\n", fold );
		std::vector<Sequence*> negSequenceSet = sampleSequenceSet( fold );

		// score negative sequences
		printf( " _________________________________\n"
				"|                                 |\n"
				"|  score negative set for fold-%d  |\n"
				"|_________________________________|\n\n", fold );
		scoreSequenceSet( motif_, bg, negSequenceSet );
		delete bg;
		delete em;
	}

	printf( " __________________________________\n"
			"|                                  |\n"
			"|  calculate precision and recall  |\n"
			"|__________________________________|\n\n" );
	calculatePR();

}
// score Global::posSequenceSet using Global::posFoldIndices and folds
void FDR::scoreSequenceSet( Motif* motif, BackgroundModel* bg, int testFold ){

	int K_motif = Global::modelOrder;
	int K_bg = Global::bgModelOrder;
	int W = motif_->getW();
	int i, LW1, y, j;
	std::vector<int> testIndices = Global::posFoldIndices[testFold];

	float maxS;									// maximal logOddsScore over all positions for each sequence

	for( size_t n = 0; n < testIndices.size(); n++ ){
		LW1 = Global::posSequenceSet->getSequences()[testIndices[n]]->getL() - W + 1;
		maxS = -FLT_MAX;						// reset maxS to a minimum negative number
		for( i = 0; i < LW1; i++ ){
			for( j = 0; j < W; j++ ){
				y = Global::posSequenceSet->getSequences()[testIndices[n]]
				                 ->extractKmer( i+j, std::min( i+j, K_motif ) );
				// TODO: This needs to be changed after the vLog() function is set up:
				posS_[testIndices[n]][i] += ( logf( motif->getV()[K_motif][y][j] )
											- logf( bg->getV()[K_bg][y] ) );
			}
			if( posS_[testIndices[n]][i] > maxS ){
				maxS = posS_[testIndices[n]][i];
			}
			posS_all_.push_back( posS_[testIndices[n]][i]);
		}
		posS_max_.push_back( maxS );
	}

}

// score sequences from negative sequence set
void FDR::scoreSequenceSet( Motif* motif, BackgroundModel* bg, std::vector<Sequence*> seqSet ){

	int K_motif = Global::modelOrder;
	int K_bg = Global::bgModelOrder;
	int W = motif_->getW();
	int i, LW1, y, j;
	float maxS;											// maximal logOddsScore over all positions for each sequence

	for( size_t n = 0; n < seqSet.size(); n++ ){
		LW1 = seqSet[n]->getL() - W + 1;
		maxS = -FLT_MAX;
		for( i = 0; i < LW1; i++ ){
			for( j = 0; j < W; j++ ){
				y = seqSet[n]->extractKmer( i+j, std::min(i+j, K_motif ) );
				negS_[n][i] += ( logf( motif->getV()[K_motif][y][j] ) - logf( bg->getV()[K_bg][y] ) );
			}
			if( negS_[n][i] > maxS ){
				maxS = negS_[n][i];
			}
			negS_all_.push_back( negS_[n][i]);
		}
		negS_max_.push_back( maxS );
	}
}

// generate negative sequences based on each test set
std::vector<Sequence*> FDR::sampleSequenceSet( int fold ){

	int k, y, yk, i, L;
	int A = Alphabet::getSize();

	// calculate (k+1)-mer frequencies for the test set:
	// 1. reset freqs_[k][y] to 0
	for( k = 0; k < Global::modelOrder + 1; k++ ){
		for( y = 0; y < ipow( A, k+1 ); y++ ){
			freqs_[k][y] = 0.0f;
		}
	}
	// 2. count (k+1)-mers y
	std::vector<Sequence*> posSeqSet = Global::posSequenceSet->getSequences();
	for( size_t idx = 0; idx < Global::posFoldIndices[fold].size(); idx++ ){
		L = posSeqSet[idx]->getL();
		for( k = 0; k < Global::modelOrder + 1; k++ ){
			for( i = 0; i < L; i++ ){
				y = posSeqSet[idx]->extractKmer( i, std::min( i, k ) );
				freqs_[k][y]++;
			}
		}
	}

	// 3. normalization of counts y
	// for k >0:
	for( k = Global::modelOrder; k > 0; k-- ){
		for( y = 0; y < ipow( A, k+1 ); y++ ){
			yk = y / A;
			freqs_[k][y] /= freqs_[k-1][yk];
		}
	}
	// for k = 0:
	float sum = 0.0f;
	for( y = 0; y < A; y++ )	sum += freqs_[0][y];
	for( y = 0; y < A; y++ )	freqs_[0][y] /= sum;

	// generate sample sequence set based on the (k+1)-mer frequencies
	std::vector<Sequence*> sampleSet;
	int posN = Global::posSequenceSet->getN();
	int negN = posN * Global::mFold / Global::cvFold;	// TODO: default Global::mFold / Global::cvFold to an integer
														// but it has to be considered if the quotient is not an integer
	for( i = 0; i < negN; i++ ){
		sampleSet.push_back( sampleSequence() );
	}

	return sampleSet;

}

// generate sample sequence based on k-mer probabilities
Sequence* FDR::sampleSequence(){

	int L = Global::posSequenceSet->getMinL();  		// TODO: better to use the average length of all the sequences
	uint8_t* sequence = ( uint8_t* )calloc( L, sizeof( uint8_t ) );
	std::string header = "sample sequence";

	int i, k, yk;
	double f, random;
	uint8_t a;
	int A = Alphabet::getSize();
	int K = Global::modelOrder;

	std::vector<int> Y;
	for( k = 0; k < Global::modelOrder + 2; k++ ){
		Y.push_back( ipow( A, k ) );
	}

	// get a random number for the first nucleotide
	random = ( double )rand() / ( double )RAND_MAX;
	f = 0.0f;
	for( a = 1; a <= A; a++ ){
		f += freqs_[0][a-1];
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
			yk += ( sequence[i-k] - 1 ) * Y[k];
		}

		// assign a nucleotide based on K-mer frequency
		f = 0.0f;
		for( a = 1; a <= A; a++ ){
			f += freqs_[std::min( i, K )][yk+a-1];
			if( random <= f ){
				sequence[i] = a;
				break;
			}
			if( sequence[i] == 0 )	sequence[i] = a;	// Trick: this is to solve the numerical problem
		}
	}

	Sequence* sampleSequence = new Sequence( sequence, L, header, Y, Global::revcomp );

	// only for testing: print out generated sequence
//	for( i = 0; i < L; i++ ){
//		std::cout << unsigned( sampleSequence->getSequence()[i] );
//	}
//	std::cout << std::endl;

	return sampleSequence;
}

void FDR::calculatePR(){

	// sort log odds scores from large to small
	// for ZOOPS model:
	std::sort( posS_max_.begin(), posS_max_.end(), std::greater<float>() );
	std::sort( negS_max_.begin(), negS_max_.end(), std::greater<float>() );
	// for MOPS model:
	std::sort( posS_all_.begin(), posS_all_.end(), std::greater<float>() );
	std::sort( negS_all_.begin(), negS_all_.end(), std::greater<float>() );

	// rank and score these log odds score values
	size_t size_max = posS_max_.size() + negS_max_.size();
	size_t size_all = posS_all_.size() + negS_all_.size();
	size_t i;
	size_t i_posM = 0, i_negM = 0;				// index for arrays storing the max log odds scores
	size_t i_posA = 0, i_negA = 0;				// index for arrays storing the complete log odds scores

	// For "many occurrence per sequence" (MOPS) model
	float max_diff = 0.0f;						// maximal difference between TFP and FP
	size_t idx_max = 0;								// index when the difference between TFP and FP reaches maximum
	for( i = 0; i < size_all; i++ ){
		if( posS_all_[i_posA] > negS_all_[i_negA] || i_posA == 0 || i_negA == negS_all_.size() ){
			i_posA++;
		} else {
			i_negA++;
		}
		FP_mops_[i] = static_cast<float>( i_negA ) / static_cast<float>( Global::mFold );
		TFP_mops_[i] = static_cast<float>( i_posA );
		if( max_diff < TFP_mops_[i] - FP_mops_[i] ){
			max_diff = TFP_mops_[i] - FP_mops_[i];
		}
		if( TFP_mops_[i] - FP_mops_[i] == max_diff ){
			idx_max = i;
		}
	}
	for( i = 0; i < idx_max; i++ ){
		P_mops_.push_back(  1 - FP_mops_[i] / TFP_mops_[i] );
		R_mops_.push_back( ( TFP_mops_[i] - FP_mops_[i] ) / max_diff );
	}

	// For "zero or one occurrence per sequence" (ZOOPS) model
	int N = Global::posSequenceSet->getN();
	for( i = 0; i < size_max; i++ ){
		if( posS_max_[i_posM] > negS_max_[i_negM] || i_posM == 0 || i_negM == negS_max_.size() ){
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

	for( i = 0; i < P_zoops_.size(); i++ ){
		ofile_zoops_p << P_zoops_[i] << std::endl;
		ofile_zoops_r << R_zoops_[i] << std::endl;
	}
	for( i = 0; i < P_mops_.size(); i++ ){
		ofile_mops_p << P_mops_[i] << std::endl;
		ofile_mops_r << R_mops_[i] << std::endl;
	}

	size_t N = Global::posSequenceSet -> getN();
	for( i = 0; i < 20 * N; i++ ){
		ofile_mops_fp << FP_mops_[i] << std::endl;
		ofile_mops_tfp << TFP_mops_[i] << std::endl;
	}
}
