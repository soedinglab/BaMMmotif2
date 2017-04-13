#ifndef MOTIF_H_
#define MOTIF_H_

#include <assert.h>
#include <math.h>	// e.g. logf

#include "Global.h"
#include "../shared/BackgroundModel.h"
#include "../shared/utils.h"

class Motif {

public:

	Motif( int length );
	Motif( const Motif& other );								// copy constructor
	~Motif();

	void initFromBaMMPattern( std::string pattern );			// initialize v from IUPAC pattern (PEnG!motif)

	void initFromBindingSites( char* filename );				// initialize v from binding site file

	void initFromPWM( float** PWM, int asize, int count );		// initialize v from PWM file

	void initFromBaMM( char* filename );						// initialize v from Bayesian Markov model file and set isInitialized

	int			getC();											// get the count of motifs N
	int         getW(); 										// get motif length w
	float***    getV();											// get conditional probabilities v
	float***	getP();											// get probabilities p
	float**		getS();											// get log odds scores for the highest order K at position j
	int***		getN();											// get the counts of (k+1)-mer for all y at motif position j
	void        updateV( float*** n, float** alpha, int order );// update v for EM
	void        updateVz_n( int*** n, float** alpha, int order );	// update v for Collapsed Gibbs sampling

	void		calculateP();									// calculate probabilities p
	void		calculateS( BackgroundModel* bg );				// calculate log odds scores for the highest order K at position j
	void		calculateLinearS( BackgroundModel* bg );		// calculate log odds scores for the highest order K at position j
																// in linear space for speeding up

	void 		print();					   					// print v to console
	void 		write( int N );					   				// write v (basename.ihbcp/.ihbp). Include header with alphabetType

private:

	bool		isInitialized_ = false;		    				// assert in all public methods

	int			C_ = 0;											// count the number of binding sites
	int 		W_;					    						// motif length
	float***    v_;				                				// conditional probabilities for (k+1)-mers y at motif position j
	float*		f_bg_;											// monomer frequencies from negative set
	float***	p_;												// probabilities for (k+1)-mers y at motif position j
	float**		s_;												// log odds scores for (K+1)-mers y at motif position j
	int***		n_;												// counts of (k+1)-mer for all y at motif position j

	void 		calculateV();									// calculate v from k-mer counts n and global alphas

	std::vector<int>	Y_;										// contains 1 at position 0
																// and the number of oligomers y for increasing order k (from 0 to K_) at positions k+1
																// e.g. alphabet size_ = 4 and K_ = 2: Y_ = 1 4 16 64
};

inline int Motif::getW(){
	return W_;
}

inline float*** Motif::getV(){
	return v_;
}

// update v from fractional k-mer counts n and current alphas (e.g for EM)
inline void Motif::updateV( float*** n, float** alpha, int K ){

	assert( isInitialized_ );

	int y, j, k, y2, yk;

	// sum up the n over (k+1)-mers at different position of motif
	std::vector<float> sumN;
	sumN.resize( W_ );
	for( j = 0; j < W_; j++ ){
		for( y = 0; y < Y_[1]; y++ ){
			sumN[j] += n[0][y][j];
		}
	}

	// for k = 0, v_ = freqs:
	for( y = 0; y < Y_[1]; y++ ){
		for( j = 0; j < W_; j++ ){
			v_[0][y][j] = ( n[0][y][j] + alpha[0][j] * f_bg_[y] )
						/ ( sumN[j] + alpha[0][j] );
		}
	}

	// for k > 0:
	for( k = 1; k < K+1; k++ ){
		for( y = 0; y < Y_[k+1]; y++ ){
			y2 = y % Y_[k];									// cut off the first nucleotide in (k+1)-mer
			yk = y / Y_[1];									// cut off the last nucleotide in (k+1)-mer
			//todo: Merge first loop into second one (by allowing contexts to extend left of the motif)
			for( j = 0; j < k; j++ ){						// when j < k, i.e. p(A|CG) = p(A|C)
				v_[k][y][j] = v_[k-1][y2][j];
			}
			//todo: Vectorize in AVX2 / SSE2
			for( j = k; j < W_; j++ ){
				v_[k][y][j] = ( n[k][y][j] + alpha[k][j] * v_[k-1][y2][j] )
							/ ( n[k-1][yk][j-1] + alpha[k][j] );
			}
		}
	}
}


// update v from integral k-mer counts n and current alphas (e.g for CGS)
inline void Motif::updateVz_n( int*** n, float** alpha, int K ){

	assert( isInitialized_ );

	int y, j, k, y2, yk;

	// sum up the n over (k+1)-mers at different position of motif
	std::vector<int> sumN;
	sumN.resize( W_ );
	for( j = 0; j < W_; j++ ){
		for( y = 0; y < Y_[1]; y++ ){
			sumN[j] += n[0][y][j];
		}
	}

	// for k = 0, v_ = freqs:
	for( y = 0; y < Y_[1]; y++ ){
		for( j = 0; j < W_; j++ ){
			v_[0][y][j] = ( ( float )n[0][y][j] + alpha[0][j] * f_bg_[y] )
						/ ( ( float )sumN[j] + alpha[0][j] );
		}
	}

	// for 1 <= k <= K:
	for( k = 1; k < K+1; k++ ){
		for( y = 0; y < Y_[k+1]; y++ ){
			y2 = y % Y_[k];									// cut off the first nucleotide in (k+1)-mer
			yk = y / Y_[1];									// cut off the last nucleotide in (k+1)-mer
			//todo: Merge first loop into second one (by allowing contexts to extend left of the motif)
			for( j = 0; j < k; j++ ){						// when j < k, i.e. p(A|CG) = p(A|C)
				v_[k][y][j] = v_[k-1][y2][j];
			}
			//todo: Vectorize in AVX2 / SSE2
			for( j = k; j < W_; j++ ){
				v_[k][y][j] = ( ( float )n[k][y][j] + alpha[k][j] * v_[k-1][y2][j] )
							/ ( ( float )n[k-1][yk][j-1] + alpha[k][j] );
			}
		}
	}
}
#endif /* MOTIF_H_ */
