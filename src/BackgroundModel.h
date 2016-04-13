/*
 * BackgroundModel.h
 *
 *  Created on: Apr 1, 2016
 *      Author: wanwan
 */

#ifndef BACKGROUNDMODEL_H_
#define BACKGROUNDMODEL_H_

#include "Global.h"
#include "SequenceSet.h"

class BackgroundModel {

public:

	BackgroundModel(){
		// allocate v
	}
	~BackgroundModel();

	void init( std::vector< std::vector<int> > folds = NULL ){

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

	float** getV(){ return v_; };	// get conditional probabilities for k-mers

	void 	print();				// print background model to console
	void 	write();			    // write background model to file basename.bmm in output directory

private:

	float** v_;					    // conditional probabilities for k-mers

	void 	calculateV( int** n );  // calculate v from k-mer counts n
};


#endif /* BACKGROUNDMODEL_H_ */
