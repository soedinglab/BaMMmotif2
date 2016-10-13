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
	std::vector<int>	trainFolds_;		// fold indices for training set

	std::vector<std::vector<float>> posSetScores_;
	std::vector<std::vector<float>> negSetScores_;

	std::vector<float>	P_zoops_;			// precision for ZOOPS model
	std::vector<float>	R_zoops_;			// recall for ZOOPS model

	std::vector<float>	P_mops_;			// precision for MOPS model
	std::vector<float>	R_mops_;			// recall for MOPS model
	float*				FP_mops_;			// false positive values for MOPS model
	float*				TFP_mops_;			// true and false positive values for MOPS model

	std::vector<int>	Y_;					// contains 1 at position 0
											// and the number of oligomers y for increasing order k (from 0 to K_) at positions k+1
											// e.g. alphabet size_ = 4 and K_ = 2: Y_ = 1 4 16 64
	float**				testsetFreq_;		// kmer frequencies in the test set
	int** 				testsetCount_;		// kmer counts in the test set

							// generate background sample sequence set based on the full positive sequence set
	std::vector<Sequence*>	sampleSequenceSet();

							// generate background sample sequence set based on each test set
	std::vector<Sequence*>	sampleSequenceSet( std::vector<Sequence*> seqSet );

							// generate negative sequence based on each sequence in the test set
	Sequence* 				sampleSequence( int length, float** freq );

							// score sequences for both positive and negative sets
	std::vector<std::vector<float>>	scoreSequenceSet( Motif* motif, BackgroundModel* bg, std::vector<Sequence*> seqSet );

							// calculate precision and recall for both ZOOPS and MOPS models
	void 		   			calculatePR();

							// calculate kmer frequencies for the test set
	void					calculateFreq( std::vector<Sequence*> seqs );

};

#endif /* FDR_H_ */
