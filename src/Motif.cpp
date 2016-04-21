/*
 * Motif.cpp
 *
 *  Created on: Apr 19, 2016
 *      Author: administrator
 */
/*
#include "Motif.h"

Motif::Motif( int length ){
	// allocate memory for v_
	length_ = length;
}

Motif::Motif( const Motif& other ){ //deep copy
	length_ = other.length_;
	if( other.v_ != NULL ){
		v_ = ( float*** )malloc( ( Global::modelOrder+1 )*sizeof( float** ) );
		for( int k = 0; k <= Global::modelOrder; k++ ){
			v_[k] = ( float** )calloc( pow( Alphabet::getSize(), k+1 ), sizeof( float* ) );
			for( int y = 0; y < pow( Alphabet::getSize(), k+1 ); y++ ){
				v_[k][y] = ( float* )calloc( length_, sizeof( float ) );
				memcpy( v_[k][y], other.v_[k][y], length_ * sizeof( float ) );
			}
		}
		isInitialized_ = true;
	} else{
		v_ = NULL;
	}

}

// initialize v from IUPAC pattern (BAMM pattern)
void Motif::initFromIUPACPattern( char* pattern ){
	// calculate k-mer counts n
	// calculate v from k-mer counts n using calculateV(n)
	// set isInitialized
}

// initialize v from binding site file
void Motif::initFromBindingSites( char* filename ){
	// scan binding sites and calculate k-mer counts n
	// calculate v from k-mer counts n using calculateV(n)
	// set isInitialized
}

// initialize v from PWM file
void Motif::initFromPWM( char* filename ){
	// set higher-order conditional probabilities to PWM probabilities
	// v[k][y][j] = PWM[0][y][j]
	// set isInitialized
}

// initialize v from Bayesian Markov model file and set isInitialized
void Motif::initFromBayesianMarkovModel( char* filename ){

}

int Motif::getLength(){
	return length_;
}

float*** Motif::getV(){
	return v_;
}

void Motif::updateV( float*** n, float** alpha ){
	assert( isInitialized_ );
	// update v from fractional k-mer counts n and current alphas
}

void Motif::print(){}
void Motif::write(){}


*/
