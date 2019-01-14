#ifndef MOTIF_H_
#define MOTIF_H_

#include <assert.h>
#include <math.h>	// e.g. logf

#include "BackgroundModel.h"
#include "../Global/utils.h"

class Motif {

public:

	Motif( size_t length );

	Motif( const Motif& other );					// copy constructor
	~Motif();

	void initFromBindingSites( char* indir, size_t l_flank, size_t r_flank );

	void initFromPWM( float** PWM, size_t asize, SequenceSet* posSet, float q );

	void initFromBaMM( char* indir, size_t l_flank, size_t r_flank );

	size_t				getW(); 					// get motif length w
	size_t				getK();						// get motif model order k
    float               getQ();                     // get estimated motif fraction on the sequences q
	float**				getA();						// get motif hyperparameter alpha
	float***    		getV();						// get conditional probabilities v
    float***    		getP();						// get probabilities v
    float**				getS();						// get log odds scores for the highest order K at position j
	std::vector<size_t> getY();

	void        		updateV( float*** n, float** alpha, size_t k );

	void				calculateP();				// calculate probabilities p
	void				calculateLogS( float** Vbg, size_t K_bg );
	void				calculateLinearS( float** Vbg, size_t K_bg ); // calculate S in linear space for speeding up

	void 				print();					// print v to console
	void 				write( char* odir, std::string basename );

private:

	bool				isInitialized_ = false;		// assert in all public methods

	size_t				C_;							// count of the sequences
	size_t 				W_;							// motif length
	size_t				K_;							// motif model order
    float_t             q_;                         // estimated motif fraction on the sequences
	float**			 	A_;							// hyperparameter alphas
	float***    		v_;				        	// conditional probabilities for (k+1)-mers y at motif position j
    float**	            v_bg_;			            // conditional probabilities of background model
    size_t              k_bg_;                      // order of background model
	float***			p_;							// probabilities for (k+1)-mers y at motif position j
	float**				s_;							// log odds scores for (K+1)-mers y at motif position j
	int***			    n_;							// fractional counts of (k+1)-mer for all y at motif position j
	std::vector<size_t>	Y_;							// contains 1 at position 0
													// and the number of oligomers y for increasing order k (from 0 to K_) at positions k+1
													// e.g. alphabet size_ = 4 and K_ = 2: Y_ = 1 4 16 64

	void 				calculateV( int*** n );	// calculate v from k-mer counts n and global alphas

};

inline size_t Motif::getW(){
	return W_;
}

inline size_t Motif::getK(){
	return K_;
}

inline float Motif::getQ(){
    return q_;
}

inline float** Motif::getA(){
	return A_;
}

inline float*** Motif::getV(){
	return v_;
}

inline float*** Motif::getP(){
    return v_;
}

inline std::vector<size_t> Motif::getY(){
	return Y_;
}

// update v from fractional k-mer counts n and current alphas
inline void Motif::updateV( float*** n, float** alpha, size_t K ){

	assert( isInitialized_ );

	// sum up the n over (k+1)-mers at different position of motif
	std::vector<float> sumN( W_ );
    for( size_t j = 0; j < W_; j++ ) {
        sumN[j] = 0.f;
    }

	for( size_t y = 0; y < Y_[1]; y++ ){
        for( size_t j = 0; j < W_; j++ ) {
            sumN[j] += n[0][y][j];
        }
	}

	// for k = 0, v_ = freqs:
	for( size_t y = 0; y < Y_[1]; y++ ){
		for( size_t j = 0; j < W_; j++ ){
			v_[0][y][j] = ( n[0][y][j] + alpha[0][j] * v_bg_[0][y] )
						/ ( sumN[j] + alpha[0][j] );
            assert( v_[0][y][j] <= 1.f );
		}
	}

	// for k > 0:
	for( size_t k = 1; k < K+1; k++ ){
		for( size_t y = 0; y < Y_[k+1]; y++ ){
			size_t y2 = y % Y_[k];				// cut off the first nucleotide
			size_t yk = y / Y_[1];				// cut off the last nucleotide
			// todo: Merge first loop into second one
            // (by allowing contexts to extend left of the motif)
			for( size_t j = 0; j < k; j++ ){	// when j < k, p(A|CG) = p(A|C)
				v_[k][y][j] = v_[k-1][y2][j];
			}
			// todo: Vectorize in AVX2 / SSE2
			for( size_t j = k; j < W_; j++ ){
				v_[k][y][j] = ( n[k][y][j] + alpha[k][j]* v_[k-1][y2][j] )
							/ ( n[k-1][yk][j-1] + alpha[k][j] );
			}
		}
	}
}

#endif /* MOTIF_H_ */
