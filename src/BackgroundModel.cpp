/*
 * BackgroundModel.cpp
 *
 *  Created on: Apr 18, 2016
 *      Author: administrator
 */

#include <numeric>
#include "BackgroundModel.h"
#include "SequenceSet.h"

BackgroundModel::BackgroundModel(){
	// allocate memory for v
	v_ = ( float** )calloc( ( Global::modelOrder+1 ), sizeof( float* ) );
	for( unsigned int k = 0; k <= Global::modelOrder; k++ ){
		v_[k] = ( float* )calloc( pow( Alphabet::getSize(), k+1 ), sizeof( float ) );
	}
}

BackgroundModel::~BackgroundModel(){
	if( v_ ) free( v_ );
}

void BackgroundModel::init( std::vector< std::vector<int> > folds ){

    if( folds == NULL ){
        std::vector<int> folds( Global::nFolds );
        std::iota( folds.begin(), folds.end(), 0 );
    }

    int** n;
    // calculate k-mer counts n from SequenceSet using Global::negFoldIndices and folds
	for( unsigned int f = 0; f < folds.size(); f++ ){
		std::vector<int> fold = folds[f];
		for( unsigned int i = 0; i < Global::negFoldIndices[f].size(); i++ ){
			Sequence sequences = Global::negSequenceSet->getSequences()[Global::negFoldIndices[f][i]];
		}
	}

    // calculate v from k-mer counts n
	calculateV( n );
}

void BackgroundModel::calculateV( int** n ){
	for( unsigned int k = 0; k < Global::modelOrder + 1; k++ ){
		for( int y = 0; y < pow( Alphabet::getSize(), k + 1 ); y++ ){
			v_[k][y] = NULL;// To be fulfilled
		}
	}
}

float** BackgroundModel::getV(){
    return v_;
}

void BackgroundModel::print(){}

void BackgroundModel::write(){}

