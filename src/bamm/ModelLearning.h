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

	ModelLearning( Motif* motif, BackgroundModel* bg, std::vector<int> folds = std::vector<int>() );
	~ModelLearning();

	int						EM();
	void 					GibbsSampling();

	void					print();
	void					write( int n );
	Motif*					getMotif();


private:

	Motif* 					motif_;				// motif to optimize within the EM
	BackgroundModel*		bg_;				// background model

	std::vector<int>		folds_;				// folds to iterate over, for cross-validation

	float**					s_;					// log scores of each (K+1)-mer at each position
	float** 				r_;		        	// responsibilities at position i in sequence n
	float*** 				n_;	            	// fractional counts n for (k+1)-mers y at motif position j
	int***					n_z_;				// n^z_j(y_1:k), the k-mer counts(for 0<k<K+2 ) with y_k's rightmost nucleotide at position j

	int*					z_;					// observed position of motif in each sequence
	float**					pos_;				// positional prior, pos[i]=0 means no motif is found on the sequence
	float** 				alpha_;	        	// pseudo-count hyper-parameter for order k and motif position j
	float 					q_; 				// hyper-parameter q specifies the fraction of sequences containing motif

	float 					llikelihood_ = 0.0f;// log likelihood for each iteration
	float					m1_t_;				// first moment for alpha optimizer
	float					m2_t_;				// second moment for alpha optimizer
	std::vector<Sequence*>	posSeqs_;			// copy positive sequences due to folds

	void 					EM_EStep();			// E-step
	void 					EM_MStep();			// M-step

	void 					EM_optimize_q();	// optimize hyper-parameter q

	void					CGS_sampling_z_q();													// sample z and q by collapsed Gibbs sampling
	inline void				CGS_updateAlphas( int order, int width, float learningrate, int t );// update alphas for all the orders up to K, given the learning rate
	float					calc_logPosterior_alphas( float** alphas, int order );				// calculate the log posterior of alphas
	float					calc_gradient_alphas( float** alphas, int order, int position );	// calculate the gradient of the log posterior of alphas

	void					testAlphaUpdate( float** alphas, int order, int width );			// only for testing, will be removed afterwards

	// test on double precision
	double					calcLogPostAlphasD( double** alphas, int order );					// calculate the log posterior of alphas
	double					calcGradLogPostAlphasD( double** alphas, int order, int position );	// calculate the gradient of the log posterior of alphas
	void					testAlphaUpdateD( double** alphas, int order, int width );			// only for testing, will be removed afterwards

	// test on the alpha optimization
	void 					debug_first_term_of_derivative();
	void 					debug_second_term_plus_prior_of_derivative();

	std::vector<int>		Y_;					// contains 1 at position 0
												// and the number of oligomers y for increasing order k (from 0 to K_) at positions k+1
												// e.g. alphabet size_ = 4 and K_ = 2: Y_ = 1 4 16 64

};


#endif /* MODELLEARNING_H_ */
