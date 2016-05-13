/*
 * EM.cpp
 *
 *  Created on: Apr 18, 2016
 *      Author: administrator
 */

#include "EM.h"

EM::EM( Motif* motif, BackgroundModel* bg ){							// for full posSet

	motif_ = motif;
	bg_ = bg;
	folds_ = std::vector<int> ( Global::nFolds );

    // allocate memory for r
	r_ = ( float** )calloc( Global::posSequenceSet->getN(), sizeof( float* ) );
    for( unsigned int n = 0; n < Global::posSequenceSet->getN(); n++ )
    	r_[n] = ( float* )calloc( Global::posSequenceSet->getSequences()[n].getL(), sizeof( float ) );

    // allocate memory for n
	n_ = ( float*** )calloc( ( Global::modelOrder+1 ), sizeof( float** ) );
	for( int k = 0; k <= Global::modelOrder; k++ ){
		n_[k] = ( float** )calloc( int( pow( Alphabet::getSize(), k+1 ) ), sizeof( float* ) );
		for( int y = 0; y < pow( Alphabet::getSize(), k+1 ); y++ )
			n_[k][y] = ( float* )calloc( motif_->getLength(), sizeof( float ) );
	}

	// allocate memory for alpha
	alpha_ = ( float** )calloc( ( Global::modelOrder+1 ), sizeof( float* ) );
		for( int k = 0; k <= Global::modelOrder; k++ )
			alpha_[k] = ( float* )calloc( motif_->getLength(), sizeof( float ) );

    // set initial values for alpha (alpha_k = beta x gamma^(order-1) with beta=20 and gamma=3) and q (0.9)
    q_ = 0.9;
    likelihood_ = 0;		// default value?
    EMIterations_ = 0;		// default value?
}

EM::EM( Motif* motif, BackgroundModel* bg, std::vector<int> folds ){	// for cross-validation
	motif_ = motif;
	bg_ = bg;
	folds_ = folds;

    // allocate memory for r
	r_ = ( float** )calloc( Global::posSequenceSet->getN(), sizeof( float* ) );
    for( unsigned int n = 0; n < Global::posSequenceSet->getN(); n++ )
    	r_[n] = ( float* )calloc( Global::posSequenceSet->getSequences()[n].getL(), sizeof( float ) );

    // allocate memory for n
	n_ = ( float*** )calloc( ( Global::modelOrder+1 ), sizeof( float** ) );
	for( int k = 0; k <= Global::modelOrder; k++ ){
		n_[k] = ( float** )calloc( int( pow( Alphabet::getSize(), k+1 ) ), sizeof( float* ) );
		for( int y = 0; y < pow( Alphabet::getSize(), k+1 ); y++ )
			n_[k][y] = ( float* )calloc( motif_->getLength(), sizeof( float ) );
	}
	// allocate memory for alpha
	alpha_ = ( float** )calloc( ( Global::modelOrder+1 ), sizeof( float* ) );
		for( int k = 0; k <= Global::modelOrder; k++ )
			alpha_[k] = ( float* )calloc( motif_->getLength(), sizeof( float ) );

    // set initial values for alpha (alpha_k = beta x gamma^(order-1) with beta=20 and gamma=3) and q (0.9)
    q_ = 0.9;
    likelihood_ = 0;		// default value?
    EMIterations_ = 0;		// default value?
}

EM::~EM(){
	free( r_ );
	free( n_ );
	free( alpha_ );
}

int EM::learnMotif(){
    // iterate over
	EStep();

	MStep();

    // * optional: optimizeAlphas()
	optimizeAlphas();

    // * optional: optimizeQ()
	optimizeQ();

    // * check likelihood for convergence
    // * remark: iterate over sequences using Global::posFoldIndices and folds_
	for( unsigned  int f=0; f < folds_.size() ; f++ ){
		int fold = folds_[f];
		for( unsigned int n=0; n < Global::posFoldIndices[fold].size(); n++ ){
			Sequence sequence = Global::posSequenceSet->getSequences()[Global::posFoldIndices[fold][n]];
		}
	}

	if( Global::verbose ) print();

    return 0;
}
void EM::EStep(){

}

void EM::MStep(){
// resetN()
// M step
}

void EM::optimizeAlphas(){
	// optimize alphas
	// motif.updateV()
}

void EM::optimizeQ(){

}

void EM::print(){

}

void EM::write(){

}
