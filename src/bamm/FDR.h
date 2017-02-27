#ifndef FDR_H_
#define FDR_H_

#include "../shared/BackgroundModel.h"
#include "../shared/utils.h"
#include "Motif.h"
#include "ModelLearning.h"
#include "SeqGenerator.h"
#include "ScoreSeqSet.h"

class FDR {

	/*
	 * FDR class deals with N-fold cross-validation:
	 * 1. calculate log odds scores for each position on each sequence;
	 * 2. generate negative sequence set based on k-mer frequencies if
	 *    no negative sequence set is given;
	 * 3. calculate true positives(TP), false positives(FP), false
	 *    discovery rate(FDR) and recall values;
	 * 4. calculate p-values due to log odds scores, in order to use
	 *    fdrtool R package for further estimation.
	 */

public:

	FDR( Motif* motif );
	~FDR();

	void 	evaluateMotif( int n );

	float 	getPrec_middle_ZOOPS();			// get precision when recall = 0.5 for ZOOPS model
	float 	getPrec_middle_MOPS();			// get precision when recall = 0.5 for MOPS model
	void 	print();
	void 	writePR( int n );
	void	writePvalues( int n );
	void 	writeLogOdds( int n );			// print out log odds scores for both positive and negative set
											// for both MOPS and ZOOPS models

											// generate a sample sequence set for simulation
											// can be separated in a different class

private:

	Motif*				motif_;				// initial motif

	std::vector<float> 	posScoreAll_;		// store log odds scores over all positions on the sequences
	std::vector<float> 	posScoreMax_;		// store maximal log odds score from each sequence
	std::vector<float> 	negScoreAll_;
	std::vector<float> 	negScoreMax_;

	std::vector<float>	ZOOPS_FDR_;			// precision for ZOOPS model
	std::vector<float>	ZOOPS_Rec_;			// recall for ZOOPS model
	std::vector<float>  ZOOPS_TP_;			// true positives for ZOOPS model
	std::vector<float>  ZOOPS_FP_;			// false positives for ZOOPS model


	std::vector<float>	MOPS_FDR_;			// precision for MOPS model
	std::vector<float>	MOPS_Rec_;			// recall for MOPS model
	std::vector<float>  MOPS_TP_;			// true positives for MOPS model
	std::vector<float>  MOPS_FP_;			// false positives for MOPS model

	float				occ_frac_;			// the fraction of motif occurrence
	float				occ_mult_;			// the number of motif occurrences per sequence

	std::vector<float>	PN_Pvalue_;			// p-values for scores from both positive and negative set with ZOOPS model
											// used for benchmarking
	std::vector<float>	ZOOPS_Pvalue_;		// p-values for scores from positive set with ZOOPS model
	std::vector<float>	MOPS_Pvalue_;		// p-values for scores from positive set with MOPS model

	std::vector<int>	Y_;					// contains 1 at position 0
											// and the number of oligomers y for increasing order k (from 0 to K_) at positions k+1
											// e.g. alphabet size_ = 4 and K_ = 2: Y_ = 1 4 16 64

							// score sequences for both positive and negative sets
	std::vector<std::vector<float>>	scoreSequenceSet( Motif* motif, BackgroundModel* bg, const std::vector<std::unique_ptr<Sequence>> & seqSet );

							// score sequences for both positive and negative sets
	std::vector<std::vector<float>>	scoreSequenceSet( Motif* motif, BackgroundModel* bg, std::vector<Sequence*> seqSet );

							// calculate precision and recall for both ZOOPS and MOPS models
	void 		   			calculatePR();

							// calculate P-values for log odds scores of positive sequences
	void					calculatePvalues();

};

#endif /* FDR_H_ */
