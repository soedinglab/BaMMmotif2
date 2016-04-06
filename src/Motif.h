/*
 * Motif.h
 *
 *  Created on: Apr 1, 2016
 *      Author: wanwan
 */

#ifndef MOTIF_H_
#define MOTIF_H_

#include "Global.h"

#include <assert.h>

class Motif {

public:

	Motif( int length ){
		// allocate memory for v_
	}
	Motif( const Motif& other ){ //deep copy
		length_ = other.length_;
		if( other.v_ != NULL ){
			v_ = ( float*** )malloc( ( Global::modelOrder+1 )*sizeof( float** ) );
			for( int k = 0; k <= Global::modelOrder; k++ ){
				v_[k] = ( float** )malloc( pow( Alphabet::getSize(), k+1 ), sizeof( float* ) );
				for( int y = 0; y < pow( Alphabet::getSize(), k+1 ) ); y++ ){
					v_[k][y] = ( float* )calloc( length_, sizeof( float ) );
					memcpy( v_[k][y], other.v_[k][y], length_ * sizeof( float ) );
				}
			}
			isInitialized_ = true;
		} else{
			v_ = NULL;
		}

	}
	~Motif();

	// initialize v from IUPAC pattern
	void initFromIUPACPattern( char* pattern ){
		// calculate k-mer counts n
		// calculate v from k-mer counts n using calculateV(n)
		// set isInitialized
	}

	// initialize v from binding site file
	void initFromBindingSites( char* filename ){
		// scan binding sites and calculate k-mer counts n
		// calculate v from k-mer counts n using calculateV(n)
		// set isInitialized
	}

	// initialize v from PWM file
	void initFromPWM( char* filename ){
		// set higher-order conditional probabilities to PWM probabilities
		// v[k][y][j] = PWM[0][y][j]
		// set isInitialized
	}

	// initialize v from inhomogeneous Markov model file and set isInitialized
	void initFromInhomogeneousMarkovModel( char* filename );

	int 		getLength(); 				// get length
	float*** 	getV();						// get v

	void updateV( float*** n, float** alpha ){
		assert( isInitialized_ );
		// update v from fractional k-mer counts n and current alphas
	}

	void 		print();					// print v to console
	void 		write();					// write v (basename.iimm)

private:

	bool		isInitialized_ = false;		// to assert in all public methods

	int 		length_;						// length of motif
	float*** 	v_[k][y][j];					// conditional probabilities for k-mers y at motif position j

	void 		calculateV( int*** n );		//calculate v from k-mer counts n and global alphas
};

#endif /* MOTIF_H_ */
