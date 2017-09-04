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
	 * and output these sequences when the log odds score
	 * is larger than certain cutoff (default:0)
	 */

public:

	ScoreSeqSet( Motif* motif,
					BackgroundModel* bg,
					std::vector<Sequence*> seqSet );
	~ScoreSeqSet();

	void score();

	std::vector<std::vector<float>> getMopsScores();
	std::vector<float> 				getZoopsScores();

	void write( char* odir, std::string basename, float cutoff, bool ss );


private:

	Motif* 							motif_;
	BackgroundModel* 				bg_;
	std::vector<Sequence*>			seqSet_;

	std::vector<std::vector<float>>	mops_scores_;
	std::vector<float>				zoops_scores_;

	std::vector<size_t>				Y_;
};


#endif /* SCORESEQSET_H_ */
