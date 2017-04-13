/*
 * ScoreSeqSet.h
 *
 *  Created on: Dec 13, 2016
 *      Author: wanwan
 */

#ifndef SCORESEQSET_H_
#define SCORESEQSET_H_

#include "Motif.h"
#include "../shared/BackgroundModel.h"

/*
 * score sequences from sequence sets
 * using learned model and background
 */

class ScoreSeqSet{

public:

	ScoreSeqSet( Motif* motif, BackgroundModel* bg, std::vector<Sequence*> seqSet );
	~ScoreSeqSet();

	void score();
	void calcPvalues( std::vector<float> neg_scores );

	std::vector<std::vector<float>> getMopsScores();
	std::vector<float> 				getScoreAll();
	std::vector<float> 				getZoopsScores();

	void write( int N, float cutoff );
	void writePvalues( int N, float cutoff );



private:

	Motif* 							motif_;
	BackgroundModel* 				bg_;
	std::vector<Sequence*>			seqSet_;

	std::vector<std::vector<float>>	mops_scores_;
	std::vector<std::vector<float>> mops_p_values_;
	std::vector<std::vector<float>> mops_e_values_;
	std::vector<float>				zoops_scores_;

	std::vector<float> 	            ScoreAll_;		 // store log odds scores over all positions on the sequences


	std::vector<int>				Y_;		         // contains 1 at position 0
											         // and the number of oligomers y for increasing order k (from 0 to K_) at positions k+1
											         // e.g. alphabet size_ = 4 and K_ = 2: Y_ = 1 4 16 64
};


#endif /* SCORESEQSET_H_ */
