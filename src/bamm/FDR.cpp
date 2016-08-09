/*
 * FDR.cpp
 *
 *  Created on: Apr 19, 2016
 *      Author: administrator
 */

#include "FDR.h"

FDR::FDR( const Motif& motif ){

	motif_ = motif;
	W_ = motif.W_;

	int posN = Global::posSequenceSet->getN();
	int negN = posN * Global::mFold / Global::cvFold;		// default Global::mFold / Global::cvFold is an integer
	int totalN = ( Global::mFold + 1 ) * Global::posSequenceSet->getN();
	int n;
	int LW1 = Global::posSequenceSet->getMaxL() - W_ + 1;

	generateFoldIndices( posN, Global::cvFold );

	// allocate memory for posLogOddsScore
	posS_all_ = ( float* )calloc( posN * LW1, sizeof( float ) );
	posS_max_ = ( float* )calloc( posN, sizeof( float* ) );

	// allocate memory for negLogOddsScore
	negS_all_ = ( float* )calloc( negN * LW1, sizeof( float ) );
	negS_max_ = ( float* )calloc( negN, sizeof( float* ) );

	// allocate memory for Precision and Recall
	P_ZOOPS_ = ( float* )calloc( posN, sizeof( float ) );
	R_ZOOPS_ = ( float* )calloc( posN, sizeof( float ) );
}

FDR::~FDR(){
	free( posS_all_ );
	free( negS_all_ );
	free( P_ZOOPS_ );
	free( R_ZOOPS_ );
}

void FDR::evaluateMotif(){

	for( int fold = 0; fold < Global::cvFold; fold++ ){

		Motif* motif = new Motif( motif_ );

		BackgroundModel* bg = new BackgroundModel( trainFoldIndices_[fold] );

		EM* em = new EM( *motif, bg, trainFoldIndices_[fold] );
		em->learnMotif();

		// score positive test sequences
		scoreSequenceSet( motif, bg, testFoldIndices_[fold] );

		// generate negative test sequences
		SequenceSet* negSequenceSet = sampleSequenceSet( em );

		// score negative test sequences
		scoreSequenceSet( motif, bg, negSequenceSet );
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
// generate fold indices for test and training sets
void FDR::generateFoldIndices( int posN, int fold ){

	testFoldIndices_.resize( fold );
	trainFoldIndices_.resize( fold );

	for( int f = 0; f < fold; f++ )
		for( int n = 0; n < posN; n++ )
			if( n % fold == f )
				testFoldIndices_[f].push_back( n );
			else
				trainFoldIndices_[f].push_back( n );
}

// score Global::posSequenceSet using Global::posFoldIndices and folds
void FDR::scoreSequenceSet( Motif* motif, BackgroundModel* bg, std::vector<int> foldIndices ){

	int n, i, LW1, k, K, y, j;
	k = Global::modelOrder;
	K = Global::bgModelOrder;
	float maxS;

	for( n = 0; n < foldIndices.size(); n++ ){
		LW1 = posSeqs_[foldIndices[n]].getL() - W_ + 1;
		maxS = 0.0f;
		for( i = 0; i < LW1; i++ ){
			for( j = 0; j < W_; j++ ){
				y = posSeqs_[foldIndices[n]].extractKmer( i + k, k );
				posS_all_[foldIndices[n]][i] += ( logf( motif->getV()[k][y][j] ) - logf( bg->getVbg()[K][y] ) );
			}
			if( posS_all_[foldIndices[n]][i] > maxS ){
				maxS = posS_all_[foldIndices[n]][i];
			}
		}
		posS_max_[foldIndices[n]] = maxS;
	}
}

// score SequenceSet sequences
void FDR::scoreSequenceSet( Motif* motif, BackgroundModel* bg, SequenceSet* seqSet ){

	int n, i, LW1, k, K, y, j;
	k = Global::modelOrder;
	K = Global::bgModelOrder;
	float maxS;

	for( n = 0; n < seqSet->getN(); n++ ){
		LW1 = seqSet->getSequences()[n].getL() - W_ + 1;
		maxS = 0.0f;
		for( i = 0; i < LW1; i++ ){
			for( j = 0; j < W_; j++ ){
				y = seqSet->getSequences()[n].extractKmer( i + k, k );
				negS_all_[n][i] += ( logf( motif->getV()[k][y][j] ) - logf( bg->getVbg()[K][y] ) );
			}
			if( negS_all_[n][i] > maxS ){
				maxS = negS_all_[n][i];
			}
		}
		negS_max_[n] = maxS;
	}
}

void FDR::calculatePR(){

	// sort log odds scores from large to small
	quickSort( posS_max_, 0, sizeof( posS_max_ ) - 1 );
	quickSort( negS_max_, 0, sizeof( negS_max_ ) - 1 );
	quickSort( posS_all_, 0, sizeof( posS_all_ ) - 1 );
	quickSort( negS_all_, 0, sizeof( negS_all_ ) - 1 );

	// rank and score these log odds score values
	int i, i_posM = 0, i_negM = 0, i_posA = 0, i_negA = 0;
	for( int i = 0; i < sizeof( posS_max_ ) + sizeof( negS_max_ ) - 1; i++ ){
		if( posS_max_[i] > negS_max_[i] || i_posM == 0 || i_negM == sizeof( negS_max_ ) ){
			i_posM ++;
		} else{
			i_negM ++;
		}
	}

}
void FDR::print(){

}
void FDR::write(){

}

