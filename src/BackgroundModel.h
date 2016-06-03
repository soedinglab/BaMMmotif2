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

	BackgroundModel();
	~BackgroundModel();

	void 	init( std::vector<int> folds );
	float* 	getV();	                		// get conditional probabilities for k-mers

	void 	print();						// print background model to console
	void 	write();			   			// write background model to file basename.bmm in output directory

private:

	void 	calculateV();  					// calculate v from k-mer counts n

	float* 	v_;					    		// conditional probabilities for k-mers
	int*	n_occ_;							// k-mer counts
	int*	powA_;							// alphabetSize to the power k
	int		Y_;								// total number of (k+1)-mers
};


#endif /* BACKGROUNDMODEL_H_ */
