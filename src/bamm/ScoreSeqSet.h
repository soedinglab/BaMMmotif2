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

	float** score();

	void write();


private:

	Motif* 					motif_;
	BackgroundModel* 		bg_;
	std::vector<Sequence*>	seqSet_;
	float**					scores_;
	std::vector<int>		Y_;					// contains 1 at position 0
												// and the number of oligomers y for increasing order k (from 0 to K_) at positions k+1
												// e.g. alphabet size_ = 4 and K_ = 2: Y_ = 1 4 16 64
};


#endif /* SCORESEQSET_H_ */
