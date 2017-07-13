/*
 * ScoreSeqSet.h
 *
 *  Created on: Dec 13, 2016
 *      Author: wanwan
 */

#ifndef SCORESEQSET_H_
#define SCORESEQSET_H_

#include "Motif.h"
#include "BackgroundModel.h"

class ScoreSeqSet{
	/*
	 * This class is aimed for:
	 * scoring sequences from the given sequence set
	 * using learned model and background model
	 */

public:

	ScoreSeqSet( Motif* motif, BackgroundModel* bg, std::vector<Sequence*> seqSet );
	~ScoreSeqSet();

	void score();

	std::vector<std::vector<float>> getMopsScores();
	std::vector<float> 				getZoopsScores();

	void write( char* dir, size_t num, float cutoff );


private:

	Motif* 							motif_;
	BackgroundModel* 				bg_;
	std::vector<Sequence*>			seqSet_;

	std::vector<std::vector<float>>	mops_scores_;
	std::vector<float>				zoops_scores_;

	std::vector<size_t>				Y_;		// contains 1 at position 0
											// and the number of oligomers y for increasing order k (from 0 to K_) at positions k+1
											// e.g. alphabet size_ = 4 and K_ = 2: Y_ = 1 4 16 64
};


#endif /* SCORESEQSET_H_ */
