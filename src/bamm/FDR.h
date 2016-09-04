#ifndef FDR_H_
#define FDR_H_

#include "../shared/BackgroundModel.h"
#include "../shared/utils.h"
#include "Motif.h"
#include "EM.h"

class FDR {

public:

	FDR( Motif* motif );
//	FDR( const Motif& motif );
	~FDR();

	void evaluateMotif();
	void print();
	void write();

private:

	Motif*				motif_;				// motif learned on full SequenceSet

	int					testFold_;			// fold index for test set
	std::vector<int>	trainFolds_;		// fold indices for training set

	float**				freqs_;				// store the frequencies of (k+1)-mer for each test set

	float**				posS_;				// store the full set of log odds scores for positive set
	float**				negS_;				// store the full set of log odds scores for negative set

	std::vector<float>	posS_all_;			// log odds scores on positive test SequenceSets
	std::vector<float> 	negS_all_;			// log odds scores on negative test SequenceSets
	std::vector<float>	posS_max_;			// largest log odds scores on positive test SequenceSet
	std::vector<float>	negS_max_;			// largest log odds scores on negative test SequenceSet

	std::vector<float>	P_zoops_;			// precision for ZOOPS model
	std::vector<float>	R_zoops_;			// recall for ZOOPS model

	std::vector<float>	P_mops_;			// precision for MOPS model
	std::vector<float>	R_mops_;			// recall for MOPS model
	float*				FP_mops_;			// false positive values for MOPS model
	float*				TFP_mops_;			// true and false positive values for MOPS model

							// generate background sample sequence set
	std::vector<Sequence*>	sampleSequenceSet( int fold );
							// generate sample sequence based on k-mer probabilities
	Sequence*				sampleSequence();

							// score Global::posSequenceSet using Global::foldIndices and folds
	void 	        		scoreSequenceSet( Motif* motif, BackgroundModel* bg, int testFold );

							// score SequenceSet sequences
	void 		    		scoreSequenceSet( Motif* motif, BackgroundModel* bg, std::vector<Sequence*> seqSet );

	void 		   			calculatePR();

};

#endif /* FDR_H_ */

