/*
 * EM.cpp
 *
 *  Created on: Apr 18, 2016
 *      Author: administrator
 */

#include <EM.h>

EM::EM( Motif* motif, BackgroundModel* bg, std::vector<int> folds = NULL ){

    if( folds == NULL ){
        std::vector<int> folds( Global::nFolds );
        std::iota( std::begin( folds ), std::end( folds ), 0 );
    }
    this.folds = folds;

    // allocate memory for r, n, and alpha
    // set initial values for alpha (alpha_k = beta x gamma^(order-1) with beta=20 and gamma=3) and q (0.9)
}

int EM::learnMotif(){
    // iterate over
    // * Estep()
    // * MStep()
    // * optional: optimizeAlphas()
    // * optional: optimizeQ()
    // * check likelihood for convergence
    // * remark: iterate over sequences using Global::posFoldIndices and folds
    //   for( f=0; f < folds.size() ; f++ ){
    //     int fold = folds[f]
    //     for( n=0; n < Global::posFoldIndices[fold].size(); n++ ){
    //       Sequence* sequence = Global::posSequenceSet.getSequences()[Global::posFoldIndices[fold][n]]
    //     }
    //   }
    // print results/statistics
    return 0;
}

void EM::print(){

}
void EM::write(){

}



