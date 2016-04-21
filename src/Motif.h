/*
 * Motif.h
 *
 *  Created on: Apr 1, 2016
 *      Author: wanwan
 */

#ifndef MOTIF_H_
#define MOTIF_H_
/*
#include <assert.h>

#include "Global.h"

class Motif {

public:

	Motif( int length );
	Motif( const Motif& other );
	~Motif();

	void initFromIUPACPattern( char* pattern );			// initialize v from IUPAC pattern  (BAMM pattern)

	void initFromBindingSites( char* filename );		// initialize v from binding site file

	void initFromPWM( char* filename );					// initialize v from PWM file

	void initFromBayesianMarkovModel( char* filename );	// initialize v from Bayesian Markov model file and set isInitialized

	int         getLength(); 							// get length
	float***    getV();									// get v

	void        updateV( float*** n, float** alpha );

	void 		print();					   			// print v to console
	void 		write();					    		// write v (basename.bmm). Include header with alphabetType

private:

	bool		isInitialized_ = false;		    		// to assert in all public methods

	int 		length_;					    		// length of motif
	float***    v_;				                		// conditional probabilities for k-mers y at motif position j, v_[k][y][j]

	void 		calculateV( int*** n );		    		// calculate v from k-mer counts n and global alphas
};
*/
#endif /* MOTIF_H_ */
