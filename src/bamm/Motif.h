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
	Motif( const Motif& other );						// copy constructor
	~Motif();

	void initFromBaMMPattern( char* pattern );			// initialize v from IUPAC pattern  (BMM pattern)

	void initFromBindingSites( char* filename );		// initialize v from binding site file

	void initFromPWM( char* filename );					// initialize v from PWM file

	void initFromBayesianMarkovModel( char* filename );	// initialize v from Bayesian Markov model file and set isInitialized

	int			getN();									// get the number of motifs N
	int         getW(); 								// get motif length w
	float***    getV();									// get conditional probabilities v
	float***	getP();									// get probabilities p
	void        updateV( float*** n, float** alpha );
    void        updateVbyK( float*** n, float** alpha, int k ); // only update V's of order k (needed for testing alpha)
	void		calculateP();							// calculate probabilities p

	void 		print();					   			// print v to console
	void 		write();					    		// write v (basename.bmm). Include header with alphabetType

private:

	bool		isInitialized_ = false;		    		// assert in all public methods

	int			N_ = 0;									// number of binding sites
	int 		W_;					    				// motif length
	float***    v_;				                		// conditional probabilities for (k+1)-mers y at motif position j
	float***	p_;										// probabilities for (k+1)-mers y at motif position j
	int***		n_;										// counts of (k+1)-mer for all y at motif position j

	void 		calculateV();							// calculate v from k-mer counts n and global alphas

	std::vector<int>	Y_;								// contains 1 at position 0
														// and the number of oligomers y for increasing order k (from 0 to K_) at positions k+1
														// e.g. alphabet size_ = 4 and K_ = 2: Y_ = 1 4 16 64
};

#endif /* MOTIF_H_ */
