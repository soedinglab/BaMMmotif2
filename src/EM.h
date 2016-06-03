/*
 * EM.h
 *
 *  Created on: Apr 1, 2016
 *      Author: wanwan
 */

#ifndef EM_H_
#define EM_H_

#include "BackgroundModel.h"
#include "MotifSet.h"

class EM {

public:

	EM( Motif* motif, BackgroundModel bg, std::vector<int> folds );
	~EM();

	int                 learnMotif();
	void                print();
	void                write();

private:

	Motif* 				motif_;			// motif to optimize within the EM
	float*				v_bg_;			// conditional probabilities for background model

	std::vector<int>	folds_;		    // folds to iterate over, for cross-validation

	float** 			r_;		        // responsibilities at position i in sequence n
	float*** 			n_;	            // fractional counts n for (k+1)-mers y at motif position j
	float** 			alpha_;	        // pseudo-count hyper-parameter for order k and motif position j
	float 				q_; 			// hyper-parameter q specifies the fraction of sequences with motif

	float 				likelihood_;	// value of Q function for current parameters
//	int 				EMIterations_;  // counter for EM iterations => can be set as a local variable
	int*				powA_;			// size of alphabet to the power k
	int					Y_;				// number of all (k+1)-mers

	void EStep();						// E-step
	void MStep();						// M-step

	void optimizeAlphas();			    // optimize alpha hyper-parameters
	void optimizeQ();					// optimize hyper-parameter q

	float sumV( Motif* motif );			// sum up the conditional probabilities
};

inline float EM::sumV( Motif* motif ){
	float sumV = 0.0f;
	for( int k = 0; k <= Global::modelOrder; k++ )
		for( int y = 0; y < powA_[k+1]; y++ )
			for( int j = 0; j < motif->getW(); j++ )
				sumV += motif->getV()[k][y][j];
	return sumV;
}

#endif /* EM_H_ */
