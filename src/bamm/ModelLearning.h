/*
 * ModelLearning.h
 *
 *  Created on: Dec 2, 2016
 *      Author: wanwan
 */

#ifndef MODELLEARNING_H_
#define MODELLEARNING_H_

#include "../shared/BackgroundModel.h"
#include "MotifSet.h"

class ModelLearning {

public:

	ModelLearning( Motif* motif, BackgroundModel* bg, std::vector<size_t> folds = std::vector<size_t>() );
	~ModelLearning();

	int						EM();
	void 					GibbsSampling();

	void					print();
	void					write( int n );
	Motif*					getMotif();


private:

	Motif* 					motif_;				// motif to optimize within the EM
	BackgroundModel*		bg_;				// background model

	std::vector<size_t>		folds_;				// folds to iterate over, for cross-validation

	float** 				r_;		        	// responsibilities at position i in sequence n
	float*** 				n_;	            	// fractional counts n for (k+1)-mers y at motif position j
	int***					n_z_;				// n^z_j(y_1:k), the k-mer counts(for 0<k<K+2 ) with y_k's
												// rightmost nucleotide at position j
	int*					z_;					// observed position of motif in each sequence
	float**					pos_;				// positional prior, pos[i]=0 means no motif is found on the sequence
	double** 				alpha_;	        	// pseudo-count hyper-parameter for order k and motif position j
	float 					q_; 				// hyper-parameter q specifies the fraction of sequences containing motif
	float 					llikelihood_ = 0.0f;// log likelihood for each iteration
	double**				m1_t_;				// first moment for alpha optimizer (ADAM)
	double**				m2_t_;				// second moment for alpha optimizer (ADAM)
	std::vector<Sequence*>	posSeqs_;			// copy positive sequences due to folds

	void 					EStep();			// E-step
	void 					MStep();			// M-step
	void 					optimize_q();		// optimize hyper-parameter q

							// sample z and q by collapsed Gibbs sampling
	void					Gibbs_sample_z_q();

							// update alphas for all the orders up to K, given the learning rate
	void					stochastic_optimize_alphas( int order, int width, float learningrate, int t );

							// calculate the gradient of the log posterior of alphas
	double					calc_gradient_alphas( double** alphas, int order, int position );

							// calculate the log posterior of a's
	double					calc_logCondProb_a( int iteration, double a, int order, int position );

							// Gibbs sampling alphas with Metropolis-Hastings algorithm
	void					GibbsMH_sample_alphas( int iteration );

							// sampling a's from the distribution of the log posterior
	void					discrete_sample_alphas( int iteration );

	// calculate the prior of alphas
	double					calc_prior_alphas( double** alphas, int order );

	// calculate the log likelihood of alphas
	double					calc_llikelihood_alphas( double** alphas, int order );

	// calculate the log posterior of alphas
	double					calc_lposterior_alphas( double** alphas, int order );

	std::vector<int>		Y_;					// contains 1 at position 0
												// and the number of oligomers y for increasing order k (from 0 to K_) at positions k+1
												// e.g. alphabet size_ = 4 and K_ = 2: Y_ = 1 4 16 64

};


#endif /* MODELLEARNING_H_ */
