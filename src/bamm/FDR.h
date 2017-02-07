#ifndef FDR_H_
#define FDR_H_

#include "../shared/BackgroundModel.h"
#include "../shared/utils.h"
#include "Motif.h"
#include "ModelLearning.h"
#include "ScoreSeqSet.h"

class FDR {

public:

	FDR( Motif* motif );
	~FDR();

	void 	evaluateMotif();

	float 	getPrec_middle_ZOOPS();			// get precision when recall = 0.5 for ZOOPS model
	float 	getPrec_middle_MOPS();			// get precision when recall = 0.5 for MOPS model
	void 	print();
	void 	write(int n);
	void 	writeLogOdds( int n);			// print out log odds scores for both positive and negative set
											// for both MOPS and ZOOPS models

private:

	Motif*				motif_;				// initial motif

	float**				testsetV_;			// k-mer frequencies in the test set
	int** 				testsetN_;			// k-mer counts in the test set

	std::vector<float> 	posScoreAll_;		// store log odds scores over all positions on the sequences
	std::vector<float> 	posScoreMax_;		// store maximal log odds score from each sequence
	std::vector<float> 	negScoreAll_;
	std::vector<float> 	negScoreMax_;

	std::vector<float>	Pre_ZOOPS_;			// precision for ZOOPS model
	std::vector<float>	Rec_ZOOPS_;			// recall for ZOOPS model
	std::vector<float>  ZOOPS_TP_;			// true positives for ZOOPS model
	std::vector<float>  ZOOPS_FP_;			// false positives for ZOOPS model
	std::vector<float>  ZOOPS_FN_;			// false negatives for ZOOPS model
	std::vector<float>  ZOOPS_TN_;			// true negatives for ZOOPS model

	std::vector<float>	Pre_MOPS_;			// precision for MOPS model
	std::vector<float>	Rec_MOPS_;			// recall for MOPS model
	std::vector<float>	FP_MOPS_;			// false positive values for MOPS model
	std::vector<float>	TFP_MOPS_;			// true and false positive values for MOPS model
	std::vector<float>  MOPS_TP_;			// true positives for MOPS model
	std::vector<float>  MOPS_FP_;			// false positives for MOPS model
	std::vector<float>  MOPS_FN_;			// false negatives for MOPS model
	std::vector<float>  MOPS_TN_;			// true negatives for MOPS model


	float 				prec_mid_ZOOPS_ = 0.0f;	// precision when recall = 0.5 for ZOOPS model
	float 				prec_mid_MOPS_ = 0.0f;	// precision when recall = 0.5 for MOPS model

	std::vector<int>	Y_;					// contains 1 at position 0
											// and the number of oligomers y for increasing order k (from 0 to K_) at positions k+1
											// e.g. alphabet size_ = 4 and K_ = 2: Y_ = 1 4 16 64

							// generate background sample sequence set based on each test set
	std::vector<std::unique_ptr<Sequence>>	sampleSequenceSet( std::vector<Sequence*> seqSet );

							// generate negative sequence based on each sequence in the test set
	std::unique_ptr<Sequence> sampleSequence( int length, float** freq );

							// score sequences for both positive and negative sets
	std::vector<std::vector<float>>	scoreSequenceSet( Motif* motif, BackgroundModel* bg, const std::vector<std::unique_ptr<Sequence>> & seqSet );

							// score sequences for both positive and negative sets
	std::vector<std::vector<float>>	scoreSequenceSet( Motif* motif, BackgroundModel* bg, std::vector<Sequence*> seqSet );

							// calculate TP FP FN and TN for ZOOPS and MOPS models
	void 					caclulateTPFPFNTN();

							// calculate precision and recall for both ZOOPS and MOPS models
	void 		   			calculatePR();

							// calculate trimer conditional probabilities for the test set
	void					calcKmerFreq( std::vector<Sequence*> seqSet );

};

#endif /* FDR_H_ */
