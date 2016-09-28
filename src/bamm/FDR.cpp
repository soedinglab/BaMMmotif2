#include "FDR.h"

#include <float.h>		// -FLT_MAX

FDR::FDR( Motif* motif ){

	motif_ = motif;

	trainFolds_.resize( Global::cvFold - 1 );
	testFold_.resize( 1 );

	for( int k = 0; k < Global::modelOrder+2; k++ ){
		Y_.push_back( ipow( Alphabet::getSize(), k ) );
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

	Motif fullMotif = Motif( *motif_ );

	for( int fold = 0; fold < Global::cvFold; fold++ ){

		// reset motif_ to the full motif
		//motif_ = Motif( fullMotif );

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

		EM* em = new EM( motif_, trainBgModel, trainFolds_ );
		em->learnMotif();

		// score positive test sequences
		printf( " _____________________________\n"
				"|                             |\n"
				"|  score test set for fold-%d  |\n"
				"|_____________________________|\n\n", fold );
		scoreSequenceSet( motif_, trainBgModel, testFold_ );

		// generate negative sequences
		printf( " __________________________________\n"
				"|                                  |\n"
				"|  create negative set for fold-%d  |\n"
				"|__________________________________|\n\n", fold );
		std::vector<Sequence*> negTestSequenceSet = sampleSequenceSet( testBgModel->getV() );
		// score negative sequences
		printf( " _________________________________\n"
				"|                                 |\n"
				"|  score negative set for fold-%d  |\n"
				"|_________________________________|\n\n", fold );
		scoreSequenceSet( motif_, trainBgModel, negTestSequenceSet );

		delete trainBgModel;
		delete em;
	}

	printf( " __________________________________\n"
			"|                                  |\n"
			"|  calculate precision and recall  |\n"
			"|__________________________________|\n\n" );
	calculatePR();

}

// score Global::posSequenceSet using Global::posFoldIndices and folds
void FDR::scoreSequenceSet( Motif* motif, BackgroundModel* bg, std::vector<int> testFold ){

	int K_motif = Global::modelOrder;
	int K_bg = Global::bgModelOrder;
	int W = motif_->getW();
	int i, LW1, k, y, j;
	std::vector<int> testIndices = Global::posFoldIndices[testFold[0]];
	float maxS;												// maximal logOddsScore over all positions for each sequence

	for( size_t n = 0; n < testIndices.size(); n++ ){
		LW1 = Global::posSequenceSet->getSequences()[testIndices[n]]->getL() - W + 1;
		maxS = -FLT_MAX;									// reset maxS to a minimum negative number
		for( i = 0; i < LW1; i++ ){
			for( j = 0; j < W; j++ ){
				k = std::min( i+j, K_motif );
				y = Global::posSequenceSet->getSequences()[testIndices[n]]->extractKmer( i+j, k );
				posScores_[testIndices[n]][i] += ( logf( motif->getV()[K_motif][y][j] )
											- logf( bg->getV()[K_bg][y] ) );
			}

			// take the complete logOddsScores for MOPS model:
			posScores_all_.push_back( posScores_[testIndices[n]][i] );

			// take the largest logOddsScore for ZOOPS model:
			if( posScores_[testIndices[n]][i] > maxS ){
				maxS = posScores_[testIndices[n]][i];
			}
		}
		posScores_max_.push_back( maxS );
	}

}

// score sequences from negative sequence set
void FDR::scoreSequenceSet( Motif* motif, BackgroundModel* bg, std::vector<Sequence*> seqSet ){

	int K_motif = Global::modelOrder;
	int K_bg = Global::bgModelOrder;
	int W = motif_->getW();
	int i, LW1, k, y, j;
	float maxS;												// maximal logOddsScore over all positions for each sequence

	for( size_t n = 0; n < seqSet.size(); n++ ){
		LW1 = seqSet[n]->getL() - W + 1;
		maxS = -FLT_MAX;
		for( i = 0; i < LW1; i++ ){
			for( j = 0; j < W; j++ ){
				k = std::min(i+j, K_motif );
				y = seqSet[n]->extractKmer( i+j, k );
				negScores_[n][i] += ( logf( motif->getV()[K_motif][y][j] ) - logf( bg->getV()[K_bg][y] ) );
			}

			// take the complete logOddsScores for MOPS model:
			negScores_all_.push_back( negScores_[n][i] );

			// take the largest logOddsScore for ZOOPS model:
			if( negScores_[n][i] > maxS ){
				maxS = negScores_[n][i];
			}
		}
		negScores_max_.push_back( maxS );
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

	int k, y, yk, i, L;
	int posN = Global::posSequenceSet->getN();
	int negN = posN * Global::mFold / Global::cvFold;	// TODO: default Global::mFold / Global::cvFold to an integer

	// generate sample sequence set based on the (k+1)-mer frequencies
	std::vector<Sequence*> sampleSet;
													// but it has to be considered if the quotient is not an integer
	for( i = 0; i < negN; i++ ){
		sampleSet.push_back( sampleSequence( v ) );
	}

	return sampleSet;

}

// generate sample sequence based on k-mer probabilities
Sequence* FDR::sampleSequence( float** v ){

	int L = Global::posSequenceSet->getMinL();  		// TODO: better to use the average length of all the sequences
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

	// only for testing: print out generated sequence
//	for( i = 0; i < L; i++ ){
//		std::cout << unsigned( sampleSequence->getSequence()[i] );
//	}
//	std::cout << std::endl;

	return sampleSequence;
}

void FDR::calculatePR(){

	// Sort log odds scores from large to small
	// for ZOOPS model:
	std::sort( posScores_max_.begin(), posScores_max_.end(), std::greater<float>() );
	std::sort( negScores_max_.begin(), negScores_max_.end(), std::greater<float>() );
	// for MOPS model:
	std::sort( posScores_all_.begin(), posScores_all_.end(), std::greater<float>() );
	std::sort( negScores_all_.begin(), negScores_all_.end(), std::greater<float>() );

	// Rank and score these log odds score values
	size_t size_max = posScores_max_.size() + negScores_max_.size();
	size_t size_all = posScores_all_.size() + negScores_all_.size();
	size_t i;
	size_t i_posM = 0, i_negM = 0;						// index for arrays storing the max log odds scores
	size_t i_posA = 0, i_negA = 0;						// index for arrays storing the complete log odds scores
	size_t N = Global::posSequenceSet->getN();

	// For "many occurrence per sequence" (MOPS) model
	float max_diff = 0.0f;								// maximal difference between TFP and FP
	size_t idx_max = 0;									// index when the difference between TFP and FP reaches maximum

	for( i = 0; i < size_all; i++ ){
		if( posScores_all_[i_posA] > negScores_all_[i_negA] ||
				i_posA == 0 || i_negA == negScores_all_.size() ){
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
		if( posScores_max_[i_posM] > negScores_max_[i_negM] ||
				i_posM == 0 || i_negM == negScores_max_.size() ){
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
