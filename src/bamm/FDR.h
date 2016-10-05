#ifndef FDR_H_
#define FDR_H_

#include "../shared/BackgroundModel.h"
#include "../shared/utils.h"
#include "Motif.h"
#include "EM.h"

class FDR {

public:

	FDR( Motif* motif );
	~FDR();

	void evaluateMotif();
	void print();
	void write();

private:

	Motif*				motif_;				// motif learned on full SequenceSet

	std::vector<int>	testFold_;			// fold index for test set
	std::vector<int>	trainFolds_;		// fold indices for training set

	float**				posScores_;			// store the full set of log odds scores for positive set
	float**				negScores_;			// store the full set of log odds scores for negative set

	std::vector<float>	posScoresAll_;		// all the log odds scores on positive test SequenceSet
	std::vector<float> 	negScoresAll_;		// all the log odds scores on negative test SequenceSet
	std::vector<float>	posScoresMax_;		// largest log odds scores on positive test SequenceSet
	std::vector<float>	negScoresMax_;		// largest log odds scores on negative test SequenceSet

	std::vector<float>	P_zoops_;			// precision for ZOOPS model
	std::vector<float>	R_zoops_;			// recall for ZOOPS model

	std::vector<float>	P_mops_;			// precision for MOPS model
	std::vector<float>	R_mops_;			// recall for MOPS model
	float*				FP_mops_;			// false positive values for MOPS model
	float*				TFP_mops_;			// true and false positive values for MOPS model

	std::vector<int>	Y_;					// contains 1 at position 0
											// and the number of oligomers y for increasing order k (from 0 to K_) at positions k+1
											// e.g. alphabet size_ = 4 and K_ = 2: Y_ = 1 4 16 64

							// generate background sample sequence set
	std::vector<Sequence*>	sampleSequenceSet( float** freq );

	Sequence* 				sampleSequence( float** freq );

							// score Global::posSequenceSet using Global::foldIndices and folds
	void 	        		scoreSequenceSet( Motif* motif, BackgroundModel* bg, std::vector<int> testFold );

							// score SequenceSet sequences
	void 		    		scoreSequenceSet( Motif* motif, BackgroundModel* bg, std::vector<Sequence*> seqSet );

							// calculate precision and recall for both ZOOPS and MOPS models
	void 		   			calculatePR();

};

#endif /* FDR_H_ */
