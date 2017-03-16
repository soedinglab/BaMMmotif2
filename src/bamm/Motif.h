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

	void initFromBaMM( char* filename );			// initialize v from Bayesian Markov model file and set isInitialized

	int			getC();											// get the count of motifs N
	int         getW(); 										// get motif length w
	float***    getV();											// get conditional probabilities v
	float***	getP();											// get probabilities p
	int***		getN();											// get the counts of (k+1)-mer for all y at motif position j
	void        updateV( float*** n, float** alpha, int order );// update v for EM
	void        updateVz_n( int*** n, float** alpha, int order );	// update v for Collapsed Gibbs sampling

	void		calculateP();									// calculate probabilities p

	void 		print();					   					// print v to console
	void 		write( int N );					   				// write v (basename.ihbcp/.ihbp). Include header with alphabetType

private:

	bool		isInitialized_ = false;		    				// assert in all public methods

	int			C_ = 0;											// count the number of binding sites
	int 		W_;					    						// motif length
	float***    v_;				                				// conditional probabilities for (k+1)-mers y at motif position j
	float***	p_;												// probabilities for (k+1)-mers y at motif position j
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

#endif /* MOTIF_H_ */
