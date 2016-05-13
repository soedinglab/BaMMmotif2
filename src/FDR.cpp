/*
 * FDR.cpp
 *
 *  Created on: Apr 19, 2016
 *      Author: administrator
 */

#include "FDR.h"

FDR::FDR( const Motif& motif ){
	motif_ = motif;
	// allocate memory for posLogOddsScore
	posS_ = ( float** )calloc( Global::posSequenceSet->getN(), sizeof( float* ) );
	for( unsigned int n = 0; n < Global::posSequenceSet->getN(); n++ )
		posS_[n] = ( float* )calloc( Global::posSequenceSet->getSequences()[n].getL(), sizeof( float ) );
	// allocate memory for negLogOddsScore
	negS_ = ( float** )calloc( Global::negSequenceSet->getN(), sizeof( float* ) );
	for( unsigned int n = 0; n < Global::negSequenceSet->getN(); n++ )
		negS_[n] = ( float* )calloc( Global::negSequenceSet->getSequences()[n].getL(), sizeof( float ) );
	// allocate memory for Precision and Recall
	P_ = ( float* )calloc( Global::posSequenceSet->getN(), sizeof( float ) );
	R_ = ( float* )calloc( Global::posSequenceSet->getN(), sizeof( float ) );
}

FDR::~FDR(){
	free( posS_ );
	free( negS_ );
	free( P_ );
	free( R_ );
}

void FDR::evaluateMotif(){

	for( unsigned int fold = 0; fold < Global::nFolds; fold++ ){

		Motif* motif = new Motif( motif_ );

		// assign training and test folds
		std::vector<int> trainingFolds;
		for( unsigned int f = 0; f < Global::nFolds; f++ ){
			if( f != fold )
				trainingFolds.push_back( f );
		}
		std::vector<int> testFold{ fold };

		BackgroundModel* bg;
		bg->init( trainingFolds );

		EM em = new EM( *motif, bg, trainingFolds );
		em.learnMotif();

		// score positive test sequences
		scoreSequenceSet( motif, bg, testFold, posS_ );

		// generate negative test sequences
		SequenceSet* negSequenceSet = sampleSequenceSet();

		// score negative test sequences
		scoreSequenceSet( motif, bg, negSequenceSet, negS_ );
	}
	calculatePR();
}

// score Global::posSequenceSet using Global::foldIndices and folds
void FDR::scoreSequenceSet( Motif* motif, BackgroundModel* bg, std::vector<int> folds, float** S ){

}

// score SequenceSet sequences
void FDR::scoreSequenceSet( Motif* motif, BackgroundModel* bg, SequenceSet* sequences, float** S ){

}

// score Sequence sequence
float* FDR::scoreSequence( Motif* motif, BackgroundModel* bg, Sequence* sequence ){
	return posS_[0]; // wrong!
}

void FDR::calculatePR(){

}
void FDR::print(){

}
void FDR::write(){

}

