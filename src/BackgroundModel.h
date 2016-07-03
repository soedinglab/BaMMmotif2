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
	float** getV();	                		// get conditional probabilities for k-mers

	void 	print();						// print background model to console
	void 	write();			   			// write background model to file basename.bmm in output directory

private:

	void 	calculateV();  					// calculate v from k-mer counts n

	int		K_;								// order of background model
	int*	powA_;							// alphabetSize to the power k
	int**	n_bg_;							// k-mer counts
	float** v_bg_;					    	// conditional probabilities for k-mers
};


#endif /* BACKGROUNDMODEL_H_ */
