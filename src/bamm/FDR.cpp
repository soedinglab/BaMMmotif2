#include "FDR.h"

FDR::FDR( Motif* motif ){

	motif_ = motif;

	trainFolds_.resize( Global::cvFold - 1 );
	testFold_.resize( 1 );

	int posN = Global::posSequenceSet->getN();
	int negN = posN * Global::mFold / Global::cvFold;		// default Global::mFold / Global::cvFold to an integer
	int LW1 = Global::posSequenceSet->getMaxL()-motif_->getW()+1;
	int n;

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

	P_ZOOPS_.reserve( posN + negN );
	R_ZOOPS_.reserve( posN + negN );
	TP_ZOOPS_ = new float[posN + negN];

	P_MOPS_.reserve( (posN + negN ) * LW1 );
	R_MOPS_.reserve( (posN + negN ) * LW1 );
	FP_MOPS_ = new float[(posN + negN ) * LW1];
	TFP_MOPS_= new float[(posN + negN ) * LW1];

}

FDR::~FDR(){
	delete[] TP_ZOOPS_;
	delete[] FP_MOPS_;
	delete[] TFP_MOPS_;
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

		BackgroundModel* bg = new BackgroundModel( Global::negSequenceSet,
													Global::bgModelOrder,
													Global::bgModelAlpha,
													Global::posFoldIndices,
													trainFolds_ );

		EM* em = new EM( motif_, bg, trainFolds_ );
		em->learnMotif();

		// score positive test sequences
		scoreSequenceSet( motif_, bg, testFold_ );

		// generate negative test sequences
		SequenceSet* negSequenceSet = sampleSequenceSet( em );

		// score negative sequences
		scoreSequenceSet( motif_, bg, negSequenceSet );

		delete bg;
	}

	calculatePR();

}

SequenceSet* FDR::sampleSequenceSet( EM* em ){

	SequenceSet* sampleSet;

	return sampleSet;

}
Sequence* FDR::sampleSequence(){

	Sequence* sequence;

	return sequence;
}

// score Global::posSequenceSet using Global::posFoldIndices and folds
void FDR::scoreSequenceSet( Motif* motif, BackgroundModel* bg, std::vector<int> testFold ){

	int i, W, LW1, k, K_bg, y, j;
	unsigned int n;
	k = Global::modelOrder;
	K_bg = Global::bgModelOrder;
	W = motif_->getW();
	float maxS;								// maximal log odds score among all the positions for each sequence
	std::vector<Sequence*> posSeqs = Global::posSequenceSet->getSequences();
	std::vector<int> testIndices = Global::posFoldIndices[testFold[0]];

	for( n = 0; n < testIndices.size(); n++ ){
		LW1 = posSeqs[testIndices[n]]->getL() - W + 1;
		maxS = 0.0f;
		for( i = 0; i < LW1; i++ ){
			for( j = 0; j < W; j++ ){
				y = posSeqs[testIndices[n]]->extractKmer( i + k, k );
				posS_[testIndices[n]][i] += ( logf( motif->getV()[k][y][j] ) - logf( bg->getVbg()[K_bg][y] ) );
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

	int n, i, W, LW1, k, K, y, j;
	k = Global::modelOrder;
	K = Global::bgModelOrder;
	W = motif_->getW();
	float maxS;

	for( n = 0; n < seqSet->getN(); n++ ){
		LW1 = seqSet->getSequences()[n]->getL() - W + 1;
		maxS = 0.0f;						// maximal log odds score among all the positions for each sequence
		for( i = 0; i < LW1; i++ ){
			for( j = 0; j < W; j++ ){
				y = seqSet->getSequences()[n]->extractKmer( i + k, k );
				negS_[n][i] += ( logf( motif->getV()[k][y][j] ) - logf( bg->getVbg()[K][y] ) );
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
	quickSort( posS_max_, 0, sizeof( posS_max_ ) - 1 );
	quickSort( negS_max_, 0, sizeof( negS_max_ ) - 1 );
	quickSort( posS_all_, 0, sizeof( posS_all_ ) - 1 );
	quickSort( negS_all_, 0, sizeof( negS_all_ ) - 1 );

	unsigned int size_max = sizeof( posS_max_ ) + sizeof( negS_max_ );
	unsigned int size_all = sizeof( posS_all_ ) + sizeof( negS_all_ );
	unsigned int i;
	unsigned int i_posM = 0, i_negM = 0;		// index for arrays storing the max log odds scores
	unsigned int i_posA = 0, i_negA = 0;		// index for arrays storing the complete log odds scores

	// rank and score these log odds score values

	// For "many occurence per sequence" (MOPS) model
	float max_diff = 0.0f;						// maximal difference between TFP and FP
	unsigned int idx_max;								// indice when the difference between TFP and FP reaches maximum
	for( i = 0; i < size_all; i++ ){
		if( posS_all_[i] > negS_all_[i] || i_posA == 0 || i_negA == sizeof( negS_all_ ) ){
			i_posA++;
		} else {
			i_negA++;
		}
		FP_MOPS_[i] = i_negA / Global::mFold;
		TFP_MOPS_[i] = i_posA;
		if( max_diff < TFP_MOPS_[i] - FP_MOPS_[i] ){
			max_diff = TFP_MOPS_[i] - FP_MOPS_[i];
		}
		if( TFP_MOPS_[i] - FP_MOPS_[i] == max_diff ){
			idx_max = i;
		}
	}
	for( i = 0; i < idx_max; i++ ){
		P_MOPS_.push_back(  1 - FP_MOPS_[i] / TFP_MOPS_[i] );
		R_MOPS_.push_back( ( TFP_MOPS_[i] - FP_MOPS_[i] ) / max_diff );
	}

	// For "zero or one occurrence per sequence" (ZOOPS) model
	int N = Global::posSequenceSet->getN();
	for( i = 0; i < size_max; i++ ){
		if( posS_max_[i] > negS_max_[i] || i_posM == 0 || i_negM == sizeof( negS_max_ ) ){
			i_posM++;
		} else {
			i_negM++;
		}
		TP_ZOOPS_[i] = i_posM;
		P_ZOOPS_.push_back( i_posM / ( i+1 ) );
		R_ZOOPS_.push_back( i_posM / N );
	}

}
void FDR::print(){

}
void FDR::write(){

}

