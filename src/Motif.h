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

	int			getN();									// get the number of motifs N
	int         getW(); 								// get motif length w
	float***    getV();									// get conditional probabilities v
	void        updateV( float*** n, float** alpha );

	void 		print();					   			// print v to console
	void 		write();					    		// write v (basename.bmm). Include header with alphabetType

private:

	bool		isInitialized_ = false;		    		// assert in all public methods
	int			N_ = 0;									// number of binding sites
	int 		w_;					    				// motif length
	float***    v_;				                		// conditional probabilities for k-mers y at motif position j
	float*		vbg_;									// conditional probabilities for monomers from background model
	int			K_ = Global::modelOrder;
	void 		calculateV( int*** n );					// calculate v from k-mer counts n and global alphas

};

#endif /* MOTIF_H_ */
