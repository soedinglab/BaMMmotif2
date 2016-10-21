#ifndef EM_H_
#define EM_H_

#include "../shared/BackgroundModel.h"
#include "MotifSet.h"

class EM {

public:

	EM( Motif* motif, BackgroundModel* bg, std::vector<int> folds = std::vector<int>() );
	~EM();

	void                    testFunctions();
	int						learnMotif();
	void					print();
	void					write();
	float 					calculateQfunc_gradient( float alpha, int order); // calculate gradient of Q-function

private:

	Motif* 					motif_;				// motif to optimize within the EM
	BackgroundModel*		bg_;				// background model
	std::vector<int>		folds_;				// folds to iterate over, for cross-validation

	float**					s_;					// log scores of each (k+1)-mer at each position
	float** 				r_;		        	// responsibilities at position i in sequence n
	float*** 				n_;	            	// fractional counts n for (k+1)-mers y at motif position j

	float** 				alpha_;	        	// pseudo-count hyper-parameter for order k and motif position j
	float 					q_ = 0.9f; 			// hyper-parameter q specifies the fraction of sequences containing motif

	float 					llikelihood_ = 0.0f;// log likelihood for each EM iteration
	unsigned int 			EMIterations_ = 0;  // counter for EM iterations
	float					Qfunc_ = 0.0f;		// Q function per each EM iteration

	std::vector<Sequence*>	posSeqs_;			// copy positive sequences due to folds

	void 					EStep();			// E-step
	void 					MStep();			// M-step

	void 					optimizeAlphas( float min = 0.1f, float max = 1e5, float tol = 0.001f);	// optimize alpha hyper-parameters
	void                    testAlphaLearning( ); // printputs and writing values to file for testing (can be deleted later)
	void 					optimizeQ();		// optimize hyper-parameter q
	float 					calculateQfunc( int k = Global::modelOrder );	// calculate incomplete Q-function
	float                   calculateLogPosterior( int k = Global::modelOrder ); // calculate log posterior likelihood
    float                   calculateLogPriors( int k = Global::modelOrder ); // calculate log prior part of log posterior



	std::vector<int>		Y_;					// contains 1 at position 0
												// and the number of oligomers y for increasing order k (from 0 to K_) at positions k+1
												// e.g. alphabet size_ = 4 and K_ = 2: Y_ = 1 4 16 64

};

#endif /* EM_H_ */
