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

	EM( Motif* motif, BackgroundModel* bg, std::vector<int> folds=NULL ){

		if( folds == NULL ){
			std::vector<int> folds( Global::nFolds );
			std::iota( std::begin( folds ), std::end( folds ), 0 );
		}
		this.folds = folds;

		// allocate memory for r, n, and alpha
		// set initial values for alpha (alpha_k = beta x gamma^(order-1) with beta=20 and gamma=3) and q (0.9)
	}
	~EM();

	int learnMotif(){
		// iterate over
		// * Estep()
		// * MStep()
		// * optional: optimizeAlphas()
		// * optional: optimizeQ()
		// * check likelihood for convergence
		// * remark: iterate over sequences using Global::posFoldIndices and folds
        //   for( f=0; f < folds.size() ; f++ ){
		//     int fold = folds[f]
		//     for( n=0; n < Global::posFoldIndices[fold].size(); n++ ){
		//       Sequence* sequence = Global::posSequenceSet.getSequences()[Global::posFoldIndices[fold][n]]
		//     }
		//   }
		// print results/statistics
	}

	void print();
	void write();

private:

	Motif* 				motif_;			// motif to optimize within the em
	BackgroundModel* 	bg_;			// background model

	std::vector<int>	folds_;			// folds to iterate over

	float** 			r_[n][i];		// responsibilities at position i in sequence n
	float*** 			n_[k][y][j];	// fractional counts n for k-mers y at motif position j

	float** 			alpha_[k][j];	// pseudocount hyperparameter for order k and motif position j
	float 				q_; 			// hyperparameter q specifies the fraction of sequences with motif

	float 				likelihood_;	// value of Q function for current parameters

	int 				EMIterations_;	// counter for EM iterations

	void EStep();						// E-step
	void MStep(){						// M-step
		// resetN()
		// M step
	}
	void optimizeAlphas(){				// optimize alpha hyper-parameters
		// optimize alphas
		// motif.updateV()
	}
	void optimizeQ();					// optimize hyper-parameter q

	void resetN();						// reset fractional counts n
};



#endif /* EM_H_ */
