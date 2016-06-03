/*
 * Motif.h
 *
 *  Created on: Apr 1, 2016
 *      Author: wanwan
 */

#ifndef MOTIF_H_
#define MOTIF_H_

#include <assert.h>

#include "Global.h"
#include "BackgroundModel.h"

class Motif {

public:

	Motif( int length );
	Motif( const Motif& other );
	~Motif();

	void initFromBaMMPattern( char* pattern );			// initialize v from IUPAC pattern  (BMM pattern)

	void initFromBindingSites( char* filename );		// initialize v from binding site file

	void initFromPWM( char* filename );					// initialize v from PWM file

	void initFromBayesianMarkovModel( char* filename );	// initialize v from Bayesian Markov model file and set isInitialized

	int         getW(); 								// get motif length w
	float***    getV();									// get conditional probabilities v
	float**		getS();									// get the log odds scores s
	void        updateV( float*** n, float** alpha );

	void 		print();					   			// print v to console
	void 		write();					    		// write v (basename.bmm). Include header with alphabetType

private:

	bool		isInitialized_ = false;		    		// to assert in all public methods

	int 		w_;					    				// motif length
	float***    v_;				                		// conditional probabilities for k-mers y at motif position j
	float**		s_;										// log odds scores of each (k+1)-mer at each position
	int*		powA_;									// size of alphabet to the power k
	int			Y_;										// number of all (k+1)-mers
	void 		calculateV_S( int*** n, int N );		// calculate v and s from k-mer counts n and global alphas, N is number of motifs

};


#endif /* MOTIF_H_ */
