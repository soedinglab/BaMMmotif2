/*
 * FDR.h
 *
 *  Created on: Apr 1, 2016
 *      Author: wanwan
 */

#ifndef FDR_H_
#define FDR_H_

#include "BackgroundModel.h"
#include "Motif.h"

class FDR {

public:

	FDR( const Motif& motif ){
		// allocate memory for posLogOddsScore, negLogOddsScore, P, and R
	}
	~FDR();

	void evaluateMotif();
	void print();
	void write();

private:

	const Motif&	motif_;				// motif learned on full SequenceSet

	float** 		posS_;		        // log-odds scores on positive test SequenceSets, posS_[n][i]
	float** 		negS_;		        // log-odds scores on negative test SequenceSets, negS_[nm][i]

	float* 			P_;					// precision
	float* 			R_;					// recall

	SequenceSet*	sampleSequenceSet();
	Sequence*		sampleSequence();

					// score Global::posSequenceSet using Global::foldIndices and folds
	void 	        scoreSequenceSet( Motif* motif, BackgroundModel* bg, std::vector<int> folds, float** S );
					// score SequenceSet sequences
	void 		    scoreSequenceSet( Motif* motif, BackgroundModel* bg, SequenceSet* sequences, float** S );
					// score Sequence sequence
	float*	        scoreSequence( Motif* motif, BackgroundModel* bg, Sequence* sequence );

	void 		    calculatePR();
};



#endif /* FDR_H_ */
