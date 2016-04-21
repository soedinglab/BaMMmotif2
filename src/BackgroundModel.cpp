/*
 * BackgroundModel.cpp
 *
 *  Created on: Apr 18, 2016
 *      Author: administrator
 */
/*
#include <numeric>
#include "BackgroundModel.h"
#include "SequenceSet.h"

BackgroundModel::BackgroundModel(){
		// allocate v
}

void BackgroundModel::init( std::vector< std::vector<int> > folds ){

    if( folds == NULL ){
        std::vector<int> folds( Global::nFolds );
        std::iota( folds.begin(), folds.end(), 0 );
    }

    // calculate k-mer counts n from SequenceSet using Global::negFoldIndices and folds
	for( unsigned int f = 0; f < folds.size(); f++ ){
		std::vector<int> fold = folds[f];
		for( unsigned int n = 0; n < Global::negFoldIndices[fold].size(); n++ ){
			Sequence* sequences = sequences_.getSequences()[Global::negFoldIndices[fold][n]];
		}
	}
    // calculate v from k-mer counts n
}

void BackgroundModel::calculateV( int** n ){

}

float** BackgroundModel::getV(){
    return v_;
}

void BackgroundModel::print(){

}

void BackgroundModel::write(){

}
*/
