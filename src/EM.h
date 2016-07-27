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

friend class Sequence;
friend class SequenceSet;

public:

	EM( Motif* motif, BackgroundModel* bg, std::vector<int> folds = std::vector<int> () );
	~EM();

	int                 learnMotif();
	void                print();
	void                write();

private:

	Motif* 					motif_;				// motif to optimize within the EM
	BackgroundModel*		bg_;				// background model

	std::vector<int>		folds_;		    	// folds to iterate over, for cross-validation

	float**					s_;					// log odds scores of each (k+1)-mer at each position
	float** 				r_;		        	// responsibilities at position i in sequence n
	float*** 				n_;	            	// fractional counts n for (k+1)-mers y at motif position j
	float***				p_;					// probabilities of PWM
	float** 				alpha_;	        	// pseudo-count hyper-parameter for order k and motif position j
	float 					q_ = 0.9f; 			// hyper-parameter q specifies the fraction of sequences containing motif

	float 					likelihood_ = 0;	// value of Q function for current parameters
	unsigned int 			EMIterations_ = 0;  // counter for EM iterations

	int						posSetN_ = Global::posSequenceSet->getN();
	std::vector<Sequence>	posSeqs_ = Global::posSequenceSet->getSequences();
	int						K_ = Global::modelOrder;
	int						k_bg_ = Global::bgModelOrder;
	int 					W_;					// motif length
	float***				v_motif_;			// conditional probabilities for motif
	float**					v_bg_;				// conditional probabilities for background model

	void EStep();								// E-step
	void MStep();								// M-step

	void optimizeAlphas();			    		// optimize alpha hyper-parameters
	void optimizeQ();							// optimize hyper-parameter q

};

#endif /* EM_H_ */
