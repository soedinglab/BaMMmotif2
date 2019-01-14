/*
 * ScoreSeqSet.h
 *
 *  Created on: Dec 13, 2016
 *      Author: wanwan
 */

#ifndef SCORESEQSET_H_
#define SCORESEQSET_H_

#include "../init/Motif.h"
#include "../init/BackgroundModel.h"

class ScoreSeqSet{
	/*
	 * This class is aimed for:
	 * scoring sequences from the given sequence set
	 * using the (learned) model and background model,
	 * find the occurrences of motifs on each sequence
	 * and output these sequences when the p-value
	 * is smaller than certain cutoff (default:0.0001)
	 */

public:

	ScoreSeqSet( Motif* motif, BackgroundModel* bg, std::vector<Sequence*> seqSet );
	~ScoreSeqSet();

	void calcLogOdds();
	void calcPvalues( std::vector<std::vector<float>> pos_mops_scores, std::vector<float> neg_all_scores );

	std::vector<std::vector<float>> getMopsScores();
	std::vector<float> 				getZoopsScores();

	void write( char* odir, std::string basename, float pvalCutoff, bool ss );
    void writeLogOdds( char* odir, std::string basename, bool ss );
    void printLogOdds();

private:

	Motif* 							motif_;
	BackgroundModel* 				bg_;
	std::vector<Sequence*>			seqSet_;

    std::vector<float>				zoops_scores_;
    std::vector<size_t>             seql_;      // create a vector to store the cumulative positive sequence lengths
	std::vector<std::vector<float>>	mops_scores_;
    std::vector<std::vector<float>> mops_p_values_;
    std::vector<std::vector<float>> mops_e_values_;
    std::vector<size_t>             z_;

    bool                            pval_is_calulated_;
	std::vector<size_t>				Y_;
};


#endif /* SCORESEQSET_H_ */
