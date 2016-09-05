#ifndef EM_H_
#define EM_H_

#include "../shared/BackgroundModel.h"
#include "MotifSet.h"

class EM {

public:

	EM( Motif* motif, BackgroundModel* bg, std::vector<int> folds = std::vector<int>() );
	~EM();

	int						learnMotif();
	void					print();
	void					write();

private:

	Motif* 					motif_;				// motif to optimize within the EM
	BackgroundModel*		bg_;				// background model
	std::vector<int>		folds_;				// folds to iterate over, for cross-validation

	float**					s_;					// log scores of each (k+1)-mer at each position
	float** 				r_;		        	// responsibilities at position i in sequence n
	float*** 				n_;	            	// fractional counts n for (k+1)-mers y at motif position j
	float***				probs_;				// probabilities of PWM

	float** 				alpha_;	        	// pseudo-count hyper-parameter for order k and motif position j
	float 					q_ = 0.9f; 			// hyper-parameter q specifies the fraction of sequences containing motif

	float 					llikelihood_ = 0.0f;// log likelihood for each EM iteration
	unsigned int 			EMIterations_ = 0;  // counter for EM iterations
	float					Qfunc_ = 0.0f;		// Q function per each EM iteration

	std::vector<Sequence*>	posSeqs_;			// copy positive sequences due to folds

	void EStep();								// E-step
	void MStep();								// M-step

	void optimizeAlphas();			    		// optimize alpha hyper-parameters
	void optimizeQ();							// optimize hyper-parameter q
	float calculateQfunc();						// calculate incomplete Q-function

};

#endif /* EM_H_ */
