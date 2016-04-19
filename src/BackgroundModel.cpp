/*
 * BackgroundModel.cpp
 *
 *  Created on: Apr 18, 2016
 *      Author: administrator
 */

#include "BackgroundModel.h"

BackgroundModel::BackgroundModel(){
		// allocate v
	}

void BackgroundModel::init( std::vector< std::vector<int> > folds = NULL ){

    if( folds == NULL ){
        std::vector<int> folds( Global::nFolds );
        std::iota( std::begin( folds ), std::end( folds ), 0 );
    }

    // calculate k-mer counts n from SequenceSet using Global::negFoldIndices and folds
    // for( f=0; f < folds.size() ; f++ ){
    //   int fold = folds[f]
    //   for( n=0; n < Global::negFoldIndices[fold].size(); n++ ){
    //     Sequence* sequence = sequences.getSequences()[Global::negFoldIndices[fold][n]]
    //   }
    // }
    // calculate v from k-mer counts n
}

float** BackgroundModel::getV(){
    return v_;
}

void BackgroundModel::print(){

}

void BackgroundModel::write(){

}
