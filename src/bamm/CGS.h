#ifndef CGS_H_
#define CGS_H_

#include "../shared/BackgroundModel.h"
#include "MotifSet.h"

class CGS {

public:
	CGS( Motif* motif, BackgroundModel* bg, std::vector<int> folds = std::vector<int>() );
	~CGS();

	void 		GibbsSampling();
	void		print();
	void		write();

private:

	Motif*				motif_;
	BackgroundModel*	bg_;

	float**				r_;					// responsibilities, posterior probability distribution
	float**				s_;					// log scores of each (K+1)-mer at each position
	int***				n_z_;				// n^z_j(y_1:k), the k-mer counts(for 0<k<K+2 ) with y_k's rightmost nucleotide at position j
	float**				alpha_;				// pseudo-count hyper-parameter for order k and motif position j
	double				q_ = 1.0f; 		// hyper-parameter q specifies the fraction of sequences containing motif
	float*				pos_;				// positional prior
	int*				z_;					// vector with coordinates z_n, observed position of motif in each sequence

	int 				CGSIterations_ = 0;

	std::vector<int>	Y_;					// contains 1 at position 0
											// and the number of oligomers y for increasing order k (from 0 to K_) at positions k+1
											// e.g. alphabet size_ = 4 and K_ = 2: Y_ = 1 4 16 64

	void		sampling_z_q( int i );		// sample z and q bt collapsed Gibbs sampling
	void		updateAlphas();				// update alphas for all the orders up to K
	float		calcGrad_logPostAlphas( float alpha, int order, int j );	// calculate the gradient of the log posterior of alphas
};

#endif  /*CGS_H_*/
