#ifndef EM_H_
#define EM_H_

#include "../shared/BackgroundModel.h"
#include "MotifSet.h"

class EM {

public:

	EM( Motif* motif, BackgroundModel* bg, std::vector<int> folds = std::vector<int>() );
	~EM();

	int                 learnMotif();
	void                print();
	void                write();

private:

	Motif* 					motif_;				// motif to optimize within the EM
	BackgroundModel*		bg_;				// background model

	std::vector<int>		folds_;				// folds to iterate over, for cross-validation

	float**					s_;					// log odds scores of each (k+1)-mer at each position
	float** 				r_;		        	// responsibilities at position i in sequence n
	float*** 				n_;	            	// fractional counts n for (k+1)-mers y at motif position j
	float***				probs_;				// probabilities of PWM
//	float**					freqs_;				// frequencies of 3-mers (for FDR sample generation)
	float** 				alpha_;	        	// pseudo-count hyper-parameter for order k and motif position j
	float 					q_ = 0.9f; 			// hyper-parameter q specifies the fraction of sequences containing motif

	float 					llikelihood_ = 0.0f;// log likelihood for each EM iteration
	unsigned int 			EMIterations_ = 0;  // counter for EM iterations
	float					Qfunc_ = 0.0f;		// Q function per each EM iteration

	unsigned int			posSetN_;
	std::vector<Sequence*>	posSeqs_;
	int						K_;
	int						k_bg_;
	int 					W_;					// motif length
	float***				v_motif_;			// conditional probabilities for motif
	float**					v_bg_;				// conditional probabilities for background model
	std::vector<int>		Y_;					// contains 1 at position 0
												// and the number of oligomers y for increasing order k (from 0 to K_) at positions k+1
												// e.g. alphabet size_ = 4 and K_ = 2: Y_ = 1 4 16 64

	void EStep();								// E-step
	void MStep();								// M-step

	void optimizeAlphas();              		// optimize alpha hyper-parameters
	void optimizeQ();							// optimize hyper-parameter q
	float calculateQfunc();						// calculate incomplete Q-function
	double calculateQfunc_gradient(double alpha ); // calculate gradient of Q-function

};

#endif /* EM_H_ */
