#include "FDR.h"

FDR::FDR( Motif* motif ){

	motif_ = motif;

	trainFolds_.resize( Global::cvFold - 1 );
	testFold_.resize( 1 );

	int posN = Global::posSequenceSet->getN();
	int negN = posN * Global::mFold / Global::cvFold;		// TODO: default Global::mFold / Global::cvFold to an integer
															// but it has to be considered if the quotient is not an integer
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

	posS_all_.reserve( posN * LW1 );
	negS_all_.reserve( negN * LW1 );
	posS_max_.reserve( posN );
	negS_max_.reserve( negN );

	P_zoops_.reserve( posN + negN );
	R_zoops_.reserve( posN + negN );

	P_mops_.reserve( ( posN + negN ) * LW1 );
	R_mops_.reserve( ( posN + negN ) * LW1 );
	FP_mops_ = new float[( posN + negN ) * LW1];
	TFP_mops_= new float[( posN + negN ) * LW1];

}

FDR::~FDR(){

	delete[] FP_mops_;
	delete[] TFP_mops_;
}

void FDR::evaluateMotif(){

	for( int fold = 0; fold < Global::cvFold; fold++ ){

		// assign training and test folds
		trainFolds_.clear();
		testFold_.clear();

		for( int f=0; f < Global::cvFold; f++ ){
			if( f != fold ){
				trainFolds_.push_back( f );
			}
		}
		testFold_.push_back( fold );

		BackgroundModel* bg = new BackgroundModel( *Global::negSequenceSet,
													Global::bgModelOrder,
													Global::bgModelAlpha,
													Global::interpolateBG,
													Global::posFoldIndices,
													trainFolds_ );

		EM* em = new EM( motif_, bg, trainFolds_ );
		em->learnMotif();

		// score positive test sequences
		scoreSequenceSet( motif_, bg, testFold_ );

		// generate negative test sequences
		SequenceSet* negSequenceSet = sampleSequenceSet( fold );

		// score negative sequences
		scoreSequenceSet( motif_, bg, negSequenceSet );

//		delete bg;
//		delete em;
	}

	calculatePR();

}

SequenceSet* FDR::sampleSequenceSet( int fold ){

	int k, y, i, L;
	// calculate (k+1)-mer frequencies for the given test fold
	std::vector<Sequence*> posSeqSet = Global::posSequenceSet->getSequences();
	for( size_t idx = 0; idx < Global::posFoldIndices[fold].size(); idx++ ){
		Sequence* posSeq = posSeqSet[idx];
		L = posSeq->getL();
		for( k = 0; k < Global::modelOrder + 1; k++ ){
			for( i = 0; i < L; i++ ){
				y = posSeq->extractKmer( i, k );
				freqs_[k][y]++;
			}
		}
	}
	for( k = 0; k < Global::modelOrder + 1; k++ ){
		float sum = 0.0f;
		for( y = 0; y < ipow(Alphabet::getSize(), k+1 ); y++ ){
			sum += freqs_[k][y];
		}
		freqs_[k][y] /= sum;
	}

	// generate sample sequence set based on the k-mer frequencies
	SequenceSet* sampleSet;

	int posN = Global::posSequenceSet->getN();
	int negN = posN * Global::mFold / Global::cvFold;
		for( i = 0; i < negN; i++ ){
			sampleSet->getSequences().push_back( sampleSequence() );
		}

	return sampleSet;

}

// generate sample sequence based on k-mer probabilities
Sequence* FDR::sampleSequence(){

	size_t L = Global::posSequenceSet->getMinL();  // TODO: better to use the average length of all the sequences
	uint8_t* sequence = ( uint8_t* )calloc( L, sizeof( uint8_t ) );
	std::string header = "sample sequence";

	int i, a, k, y;
	int A = Alphabet::getSize();
	int K = Global::modelOrder;

	// get a random number for the first nucleotide
	double random = ( double )rand() / RAND_MAX;
	float f = 0.0f;
	for( a = 0; a < A; a++ ){
		f += freqs_[0][a];
		if( random < f ){
			sequence[0] = a;
			break;
		} else {
			continue;
		}
	}

	for( i = 1; i < L; i++ ){
		// get a random double number
		random = ( double )rand() / RAND_MAX;

		// calculate y of (K+1)-mer
		y = 0;
		for( k = std::min( K, i ); k > 0; k-- ){
			y += ( sequence[i-k] - 1 ) * ipow( A, k );
		}
		y += a;

		// assign a nucleotide based on K-mer frequency
		f = 0.0f;				// reset f
		for( a = 0; a < A; a++ ){
			f += freqs_[std::min( K, i )][y];
			if( random < f ){
				sequence[i] = a;
				break;
			} else {
				continue;
			}
		}
	}

	Sequence* sampleSequence( sequence, L, header );

	return sampleSequence;
}

// score Global::posSequenceSet using Global::posFoldIndices and folds
void FDR::scoreSequenceSet( Motif* motif, BackgroundModel* bg, std::vector<int> testFold ){

	int i, W, LW1, k, K_bg, y, j;
	k = Global::modelOrder;
	K_bg = Global::bgModelOrder;
	W = motif_->getW();
	float maxS;									// maximal log odds score among all the positions for each sequence
	std::vector<Sequence*> posSeqs = Global::posSequenceSet->getSequences();
	std::vector<int> testIndices = Global::posFoldIndices[testFold[0]];

	for( size_t n = 0; n < testIndices.size(); n++ ){
		LW1 = posSeqs[testIndices[n]]->getL() - W + 1;
		maxS = 0.0f;
		for( i = 0; i < LW1; i++ ){
			for( j = 0; j < W; j++ ){
				y = posSeqs[testIndices[n]]->extractKmer( i + k, k );
				// TODO: This needs to be changed after the vLog() function is set up:
				posS_[testIndices[n]][i] += ( logf( motif->getV()[k][y][j] ) - logf( bg->getV()[K_bg][y] ) );
			}
			if( posS_[testIndices[n]][i] > maxS ){
				maxS = posS_[testIndices[n]][i];
			}
			posS_all_.push_back( posS_[testIndices[n]][i]);
		}
		posS_max_.push_back( maxS );
	}
}

// score SequenceSet sequences
void FDR::scoreSequenceSet( Motif* motif, BackgroundModel* bg, SequenceSet* seqSet ){

	int n, i, W, LW1, k, K_bg, y, j;
	k = Global::modelOrder;
	K_bg = Global::bgModelOrder;
	W = motif_->getW();
	float maxS;

	for( n = 0; n < seqSet->getN(); n++ ){
		LW1 = seqSet->getSequences()[n]->getL() - W + 1;
		maxS = 0.0f;							// maximal log odds score among all the positions for each sequence
		for( i = 0; i < LW1; i++ ){
			for( j = 0; j < W; j++ ){
				y = seqSet->getSequences()[n]->extractKmer( i + k, k );
				negS_[n][i] += ( logf( motif->getV()[k][y][j] ) - logf( bg->getV()[K_bg][y] ) );
			}
			if( negS_[n][i] > maxS ){
				maxS = negS_[n][i];
			}
			negS_all_.push_back( negS_[n][i]);
		}
		negS_max_.push_back( maxS );
	}
}

void FDR::calculatePR(){

	// sort log odds scores from large to small
	std::sort( posS_max_.begin(), posS_max_.end(), std::greater<float>() );
	std::sort( negS_max_.begin(), posS_max_.end(), std::greater<float>() );
	std::sort( posS_all_.begin(), posS_max_.end(), std::greater<float>() );
	std::sort( negS_all_.begin(), posS_max_.end(), std::greater<float>() );

	// rank and score these log odds score values
	size_t size_max = sizeof( posS_max_ ) + sizeof( negS_max_ );
	size_t size_all = sizeof( posS_all_ ) + sizeof( negS_all_ );
	size_t i;
	size_t i_posM = 0, i_negM = 0;				// index for arrays storing the max log odds scores
	size_t i_posA = 0, i_negA = 0;				// index for arrays storing the complete log odds scores

	// For "many occurrence per sequence" (MOPS) model
	float max_diff = 0.0f;						// maximal difference between TFP and FP
	size_t idx_max;								// index when the difference between TFP and FP reaches maximum
	for( i = 0; i < size_all; i++ ){
		if( posS_all_[i] > negS_all_[i] || i_posA == 0 || i_negA == sizeof( negS_all_ ) ){
			i_posA++;
		} else {
			i_negA++;
		}
		FP_mops_[i] = i_negA / Global::mFold;
		TFP_mops_[i] = i_posA;
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
		if( posS_max_[i] > negS_max_[i] || i_posM == 0 || i_negM == sizeof( negS_max_ ) ){
			i_posM++;
		} else {
			i_negM++;
		}

		P_zoops_.push_back( i_posM / ( i+1 ) );  // i_posM is the true positive value
		R_zoops_.push_back( i_posM / N );
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
	for( i = 0; i < sizeof( FP_mops_ ); i++ ){
		ofile_mops_p << FP_mops_[i] << std::endl;
		ofile_mops_r << TFP_mops_[i] << std::endl;
	}
}

