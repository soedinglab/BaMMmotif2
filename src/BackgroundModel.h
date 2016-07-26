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

	BackgroundModel( std::vector<int> folds = std::vector<int> () );
	~BackgroundModel();

	float** getVbg();	                				// get conditional probabilities for k-mers

	void 	print();									// print background model to console
	void 	write();			   						// write background model to file basename.bmm in output directory

private:

	void 	calculateVbg();  							// calculate v from k-mer counts n

	int**	n_bg_;										// k-mer counts
	float** v_bg_;					    				// conditional probabilities for k-mers
	int		K_ = Global::bgModelOrder;
	int 	negSetN_ = Global::negSequenceSet->getN();	// count of negative sequences
	std::vector<Sequence> negSeqs_ = Global::negSequenceSet->getSequences();

};

#endif /* BACKGROUNDMODEL_H_ */
