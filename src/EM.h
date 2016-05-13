/*
 * EM.h
 *
 *  Created on: Apr 1, 2016
 *      Author: wanwan
 */

#ifndef EM_H_
#define EM_H_

#include "BackgroundModel.h"
#include "Motif.h"

class EM {

public:
	EM( Motif* motif, BackgroundModel* bg );							// for full posSet
	EM( Motif* motif, BackgroundModel* bg, std::vector<int> folds );	// for cross-validation
	~EM();

	int                 learnMotif();
	void                print();
	void                write();

private:

	Motif* 				motif_;			// motif to optimize within the em
	BackgroundModel*	bg_;			// background model

	std::vector<int>	folds_;		    // folds to iterate over, for cross-validation

	float** 			r_;		        // responsibilities at position i in sequence n, r_[n][i]
	float*** 			n_;	            // fractional counts n for k-mers y at motif position j, n_[k][y][j]

	float** 			alpha_;	        // pseudocount hyperparameter for order k and motif position j, alpha_[k][j]
	float 				q_; 			// hyperparameter q specifies the fraction of sequences with motif

	float 				likelihood_;	// value of Q function for current parameters

	int 				EMIterations_;  // counter for EM iterations

	void EStep();						// E-step
	void MStep();						// M-step

	void optimizeAlphas();			    // optimize alpha hyper-parameters
	void optimizeQ();					// optimize hyper-parameter q

	void resetN();						// reset fractional counts n
};

#endif /* EM_H_ */
