/*
 * FDR.cpp
 *
 *  Created on: Apr 19, 2016
 *      Author: administrator
 */
/*
#include <FDR.h>

void FDR::evaluateMotif(){

	for( int fold = 0; fold < Global::nFolds; fold++ ){

		Motif* motif = new Motif( this.motif );

		// assign training and test folds
		std::vector<int> trainingFold;
		for( int f = 0; f < Global::nFolds; f++ ){
			if( f != fold ){
				trainingFold.push_back( f );
			}
		}
		std::vector<int> testFold{ fold };

		BackgroundModel bg;
		bg.init( trainingFold );

		EM em = new EM( *motif, bg, trainingFold );
		em.learnMotif();

		// score positive test sequences
		scoreSequenceSet( motif, bg, testFold, posS );

		// generate negative test sequences
		SequenceSet* negSequenceSet = sampleSequenceSet();
		// score negative test sequences
		scoreSequenceSet( motif, bg, negSequenceSet, negS );
	}
	calculatePR();
}

void FDR::print(){};
void FDR::write(){};

*/
