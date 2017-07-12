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
	void					write( size_t n );
	Motif*					getMotif();


private:

	Motif* 					motif_;				// motif to optimize within the EM
	BackgroundModel*		bg_;				// background model

	size_t 					K_;					// the order of the motif model
	size_t					W_;					// the width of the motif model
	std::vector<size_t>		folds_;				// folds to iterate over, for cross-validation

	float** 				r_;		        	// responsibilities at position i in sequence n
	float**					s_;					// log odds scores
	float*** 				n_;	            	// fractional counts n for (k+1)-mers y at motif position j
	size_t*					z_;					// observed position of motif in each sequence
	float**					pos_;				// positional prior, pos[i]=0 means no motif is found on the sequence
	float** 				A_;	        		// pseudo-count hyper-parameter for order k and motif position j
	float 					q_; 				// hyper-parameter q specifies the fraction of sequences containing motif
	size_t					N0_ = 0;			// count of sequences that do not contain a motif
	float 					llikelihood_ = 0.0f;// log likelihood for each iteration
	double**				m1_t_;				// first moment for alpha optimizer (ADAM)
	double**				m2_t_;				// second moment for alpha optimizer (ADAM)
	std::vector<Sequence*>	posSeqs_;			// copy positive sequences due to folds
	std::vector<size_t>		Y_;

	void 					EStep();			// E-step
	void 					MStep();			// M-step
	void 					Optimize_q();		// optimize hyper-parameter q

							// sample motif position z by collapsed Gibbs sampling
	void					Collapsed_Gibbs_sample_z();

							// sample sequence fraction q for motif by regular Gibbs sampling
	void					Gibbs_sample_q();

							// update alphas for all the orders up to K by stochastic gradient descent
	void					Optimize_alphas_by_SGD( size_t order, size_t width, float learningrate, size_t t );

							// Gibbs sampling alphas with Metropolis-Hastings algorithm
	void					GibbsMH_sample_alphas( size_t iteration );

							// sampling a's from the distribution of the log posterior
	void					Discrete_sample_alphas( size_t iteration );

							// calculate the gradient of the log posterior of alphas
	float					calc_gradient_alphas( float** alphas, size_t order, size_t position );

							// calculate the log posterior of a's
	float					calc_logCondProb_a( size_t iteration, float a, size_t order, size_t position );

							// calculate the prior of alphas
	float					calc_prior_alphas( float** alphas, size_t order );

							// calculate the log likelihood of alphas
	float					calc_llikelihood_alphas( float** alphas, size_t order );

							// calculate the log posterior of alphas
	float					calc_lposterior_alphas( float** alphas, size_t order );

};


#endif /* MODELLEARNING_H_ */
