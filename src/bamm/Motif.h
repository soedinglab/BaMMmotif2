#ifndef MOTIF_H_
#define MOTIF_H_

#include <assert.h>
#include <math.h>	// e.g. logf

#include "Global.h"
#include "../shared/BackgroundModel.h"
#include "../shared/utils.h"

class Motif {

public:

	Motif( size_t length, size_t order = Global::modelOrder, std::vector<float> alpha = Global::modelAlpha );
	Motif( const Motif& other );								// copy constructor
	~Motif();

	void initFromBaMMPattern( std::string pattern );			// initialize v from IUPAC pattern (PEnG!motif)

	void initFromBindingSites( char* filename );				// initialize v from binding site file

	void initFromPWM( float** PWM, size_t asize, size_t count );// initialize v from PWM file

	void initFromBaMM( char* filename );						// initialize v from Bayesian Markov model file and set isInitialized

	size_t		getC();											// get the count of motifs N
	size_t		getW(); 										// get motif length w
	size_t		getK();											// get motif model order k
	float***    getV();											// get conditional probabilities v
	float***	getP();											// get probabilities p
	float**		getS();											// get log odds scores for the highest order K at position j
	size_t***	getN();											// get (k+1)-mer counts for all y at motif position j
	void        updateV( float*** n, double** alpha, size_t k );// update v for EM
	void        updateVz_n( size_t*** n, double** alpha, size_t k );// update v for Collapsed Gibbs sampling

	void		calculateP();									// calculate probabilities p
	void		calculateS( float** Vbg );						// calculate log odds scores for the highest order K at position j
	void		calculateLinearS( float** Vbg );				// calculate log odds scores for the highest order K at position j
																// in linear space for speeding up

	void 		print();					   					// print v to console
	void 		write( size_t N );					  			// write v (basename.ihbcp/.ihbp). Include header with alphabetType

private:

	bool		isInitialized_ = false;		    				// assert in all public methods

	size_t		C_ = 0;											// count the number of binding sites
	size_t 		W_;					    						// motif length
	size_t		K_;												// mortif model order
	std::vector<float> A_;										// hyperparameter alphas
	float***    v_;				                				// conditional probabilities for (k+1)-mers y at motif position j
	float*		f_bg_;											// monomer frequencies from negative set
	float***	p_;												// probabilities for (k+1)-mers y at motif position j
	float**		s_;												// log odds scores for (K+1)-mers y at motif position j
	size_t***	n_;												// counts of (k+1)-mer for all y at motif position j

	void 		calculateV();									// calculate v from k-mer counts n and global alphas

	std::vector<size_t>	Y_;										// contains 1 at position 0
																// and the number of oligomers y for increasing order k (from 0 to K_) at positions k+1
																// e.g. alphabet size_ = 4 and K_ = 2: Y_ = 1 4 16 64
};

inline size_t Motif::getW(){
	return W_;
}

inline size_t Motif::getK(){
	return K_;
}

inline float*** Motif::getV(){
	return v_;
}

// update v from fractional k-mer counts n and current alphas (e.g for EM)
inline void Motif::updateV( float*** n, double** alpha, size_t K ){

	assert( isInitialized_ );

	// sum up the n over (k+1)-mers at different position of motif
	std::vector<float> sumN;
	sumN.resize( W_ );
	for( size_t j = 0; j < W_; j++ ){
		for( size_t y = 0; y < Y_[1]; y++ ){
			sumN[j] += n[0][y][j];
		}
	}

	// for k = 0, v_ = freqs:
	for( size_t y = 0; y < Y_[1]; y++ ){
		for( size_t j = 0; j < W_; j++ ){
			v_[0][y][j] = ( n[0][y][j] + static_cast<float>( alpha[0][j] ) * f_bg_[y] )
						/ ( sumN[j] + static_cast<float>( alpha[0][j] ) );
		}
	}

	// for k > 0:
	for( size_t k = 1; k < K+1; k++ ){
		for( size_t y = 0; y < Y_[k+1]; y++ ){
			size_t y2 = y % Y_[k];									// cut off the first nucleotide in (k+1)-mer
			size_t yk = y / Y_[1];									// cut off the last nucleotide in (k+1)-mer
			//todo: Merge first loop into second one (by allowing contexts to extend left of the motif)
			for( size_t j = 0; j < k; j++ ){						// when j < k, i.e. p(A|CG) = p(A|C)
				v_[k][y][j] = v_[k-1][y2][j];
			}
			//todo: Vectorize in AVX2 / SSE2
			for( size_t j = k; j < W_; j++ ){
				v_[k][y][j] = ( n[k][y][j] + static_cast<float>( alpha[k][j] ) * v_[k-1][y2][j] )
							/ ( n[k-1][yk][j-1] + static_cast<float>( alpha[k][j] ) );
			}
		}
	}
}


// update v from integral k-mer counts n and current alphas (e.g for CGS)
inline void Motif::updateVz_n( size_t*** n, double** alpha, size_t K ){

	assert( isInitialized_ );

	// sum up the monomers
	size_t sumN = 0;
	for( size_t y = 0; y < Y_[1]; y++ ){
		sumN += n[0][y][0];
	}

	// for k = 0, v_ = freqs:
	for( size_t y = 0; y < Y_[1]; y++ ){
		for( size_t j = 0; j < W_; j++ ){
			v_[0][y][j] = ( static_cast<float>( n[0][y][j] ) + static_cast<float>( alpha[0][j] ) * f_bg_[y] )
						/ ( static_cast<float>( sumN ) + static_cast<float>( alpha[0][j] ) );
		}
	}

	// for 1 <= k <= K:
	for( size_t k = 1; k < K+1; k++ ){
		for( size_t y = 0; y < Y_[k+1]; y++ ){
			size_t y2 = y % Y_[k];									// cut off the first nucleotide in (k+1)-mer
			size_t yk = y / Y_[1];									// cut off the last nucleotide in (k+1)-mer
			//todo: Vectorize in AVX2 / SSE2
			for( size_t j = 0; j < W_; j++ ){
				v_[k][y][j] = ( static_cast<float>( n[k][y][j] ) + static_cast<float>( alpha[k][j] ) * v_[k-1][y2][j] )
							/ ( static_cast<float>( n[k-1][yk][j-1] ) + static_cast<float>( alpha[k][j] ) );
			}
		}
	}
}

#endif /* MOTIF_H_ */
