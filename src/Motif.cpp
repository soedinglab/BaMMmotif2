/*
 * Motif.cpp
 *
 *  Created on: Apr 19, 2016
 *      Author: administrator
 */

#include "Motif.h"

Motif::Motif( int length ){					// allocate memory for v
	length_ = length;
	v_ = ( float*** )calloc( ( Global::modelOrder+1 ), sizeof( float** ) );
	for( unsigned int k = 0; k <= Global::modelOrder; k++ ){
		v_[k] = ( float** )calloc( pow( Alphabet::getSize(), k+1 ), sizeof( float* ) );
		for( int y = 0; y < pow( Alphabet::getSize(), k+1 ); y++ ){
			v_[k][y] = ( float* )calloc( length_, sizeof( float ) );
		}
	}
}

Motif::Motif( const Motif& other ){ 		// deep copy
	length_ = other.length_;
	if( other.v_ != NULL ){
		v_ = ( float*** )calloc( ( Global::modelOrder+1 ), sizeof( float** ) );
		for( unsigned int k = 0; k <= Global::modelOrder; k++ ){
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

Motif::~Motif(){
	if( v_ ) free( v_ );
}


// initialize v from IUPAC pattern (BMM pattern)
void Motif::initFromIUPACPattern( char* pattern ){
	// calculate k-mer counts n
	// calculate v from k-mer counts n using calculateV(n)
	// set isInitialized
}

// initialize v from binding site file
void Motif::initFromBindingSites( char const* filename ){

	std::ifstream file( filename );				// read file
	std::string seq;							// read each sequence from single line
	int length;
	int lineNum;
	int*** n;
	initN( n );

	while( getline( file, seq ).good() ){

		length = seq.length();
		lineNum ++;

		// all the binding sites should have the same length
		if( length != length_ ){
			fprintf( stderr, "Length of sequence at line %d differs.\n"
					"Binding sites should have the same length.\n", lineNum );
			exit( -1 );
		}

		// binding sites should be longer than the order of model
		if( length <= Global::modelOrder + 1 ){
			fprintf( stderr, "Length of sequences exceeds length of model.\n" );
			exit( -1 );
		}

		// scan binding sites and calculate k-mer counts n
		for( unsigned int k = 0; k < Global::modelOrder + 1; k++ ){		// k runs over all orders
			for( int j = k; j < length_; j++ ){							// j runs over all motif position
				int y = 0;
				for( int i = 0; i <= k; i++ ){							// calculate y based on (k+1)-mer bases
					y += pow( Alphabet::getSize(), i ) * ( Alphabet::getCode( seq[j-k+i] ) - 1 );
				}
				n[k][y][j]++;
			}
		}

	}

	// calculate v from k-mer counts n using calculateV(n)
	calculateV( n, length_ * lineNum );

	// set isInitialized
	isInitialized_ = true;
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

void Motif::initN( int*** n ){
	for( unsigned int k = 0; k < Global::modelOrder + 1; k++ ){
		for( int y = 0; y < pow( Alphabet::getSize(), k + 1 ); y++ ){
			for( int j = k; j < length_; j++ ){
				n[k][y][j] = 0;
			}
		}
	}
}

void Motif::calculateV( int*** n, int sum ){

	// for k = 0:
	for( int y = 0; y < pow( Alphabet::getSize(), 1 ); y++ ){
		for( int j = 0; j < length_; j++ ){
			v_[0][y][j] = n[0][y][j] / sum;
		}
	}

	// for k > 0:
	for( unsigned int k = 1; k < Global::modelOrder + 1; k++ ){
		for( int y = 0; y < pow( Alphabet::getSize(), k + 1 ); y++ ){
			for( int j = 1; j <= length_; j++ ){
				v_[k][y][j] = ( n[k][y][j] + Global::modelAlpha[k] * v_[k-1][y][j] )
						/ ( n[k-1][y][j-1] + Global::modelAlpha[k] );
			}
		}
	}
}

float*** Motif::getV(){
	return v_;
}

void Motif::updateV( float*** n, float** alpha ){
	assert( isInitialized_ );
	// update v from fractional k-mer counts n and current alphas
}

void Motif::print(){}

void Motif::write(){
	std::string opath = Global::outputDirectory  + '/' + Global::posSequenceBasename + ".probs";
	std::ofstream ofile(opath);
	for( int j = 0; j < length_; j++ ){
		for( int k = 0; k < Alphabet::getSize(); k++ ){
			for( int y = 0; y < pow( Alphabet::getSize(), k+1 ); y++ ){
				ofile << v_[k][y][j] << '\t';
			}
			ofile << std::endl;
		}
		ofile << std::endl << std::endl;
	}
}

