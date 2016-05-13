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
	for( int k = 0; k <= Global::modelOrder; k++ ){
		v_[k] = ( float** )calloc( int( pow( Alphabet::getSize(), k+1 ) ), sizeof( float* ) );
		for( int y = 0; y < pow( Alphabet::getSize(), k+1 ); y++ ){
			v_[k][y] = ( float* )calloc( length_, sizeof( float ) );
		}
	}
}

Motif::Motif( const Motif& other ){ 		// deep copy
	length_ = other.length_;
	if( other.v_ != NULL ){
		v_ = ( float*** )calloc( ( Global::modelOrder+1 ), sizeof( float** ) );
		for( int k = 0; k <= Global::modelOrder; k++ ){
			v_[k] = ( float** )calloc( int( pow( Alphabet::getSize(), k+1 ) ), sizeof( float* ) );
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
void Motif::initFromBindingSites( char* filename ){

	// initialize n
	int*** n;
	n = new int**[Global::modelOrder + 1];
	for( int k = 0; k < Global::modelOrder + 1; k++ ){
		int sizeY = int( pow( Alphabet::getSize(), k + 1 ) );
		n[k] = new int*[sizeY];
		for( int y = 0; y < pow( Alphabet::getSize(), k + 1 ); y++ ){
			n[k][y] = new int[length_];
			for( int j = 0; j < length_; j++ ){
				n[k][y][j] = 0;
			}
		}
	}

	std::ifstream file( filename );				// read file
	std::string seq;							// read each sequence from single line
	int length;
	int lineNum = 0;

	while( getline( file, seq ).good() ){
		length = seq.length();
		lineNum ++;

		if( length != length_ ){								// all the binding sites should have the same length
			fprintf( stderr, "Error: Length of sequence at line %d differs.\n"
					"Binding sites should have the same length.\n", lineNum );
			exit( -1 );
		}
		if( length <= Global::modelOrder + 1 ){					// binding sites should be longer than the order of model
			fprintf( stderr, "Error: Length of sequence at line %d exceeds length of model.\n", lineNum );
			exit( -1 );
		}
		if( length > Global::posSequenceSet->getMinL() ){		// binding sites should be shorter than the shortest posSeq
			fprintf( stderr, "Error: Length of sequence at line %d exceeds length of posSet sequence.\n", lineNum );
			exit( -1 );
		}

		// scan binding sites and calculate k-mer counts n
		for( int k = 0; k < Global::modelOrder + 1; k++ ){				// k runs over all orders
			for( int j = k; j < length_; j++ ){							// j runs over all motif position
				int y = 0;
				for( int i = 0; i <= k; i++ ){							// calculate y based on (k+1)-mer bases
					y += int( pow( Alphabet::getSize(), i ) * ( Alphabet::getCode( seq[j-k+i] ) - 1 ) );
				}
				n[k][y][j]++;
			}
		}
	}

	// calculate v from k-mer counts n using calculateV(n)
	calculateV( n, length_ * lineNum );

	if( Global::verbose ) print();

	// set isInitialized
	isInitialized_ = true;

	// free memory
//	for( int k = 0; k < Global::modelOrder + 1; k++ ){
//		for( int y = 0; y < pow( Alphabet::getSize(), k + 1 ); y++ ){
//			delete[] n[k][y];
//		}
//		delete[] n[k];
//	}
//	delete[] n;

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

void Motif::calculateV( int*** n, int sum ){

	// for k = 0:
	for( int y = 0, k = 0; y < pow( Alphabet::getSize(), k + 1 ); y++ ){
		for( int j = 0; j < length_; j++ ){
			v_[k][y][j] = n[k][y][j] / float( sum );
		}
	}

	// for k > 0:
	for( int k = 1; k < Global::modelOrder + 1; k++ ){
		for( int y = 0; y < pow( Alphabet::getSize(), k + 1 ); y++ ){
			int y2 = y / pow( Alphabet::getSize(), k );				// cut off first base in (k+1)-mer y
			int yk = y % int( pow( Alphabet::getSize(), k ) );		// cut off last base in (k+1)-mer y
			for( int j = 1; j <= length_; j++ ){
				v_[k][y][j] = ( n[k][y][j] + Global::modelAlpha[k] * v_[k-1][y2][j] )
						/ ( n[k-1][yk][j-1] + Global::modelAlpha[k] );
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

void Motif::print(){
	printf( " ___________________________\n"
			"|                           |\n"
			"| INITIALIZED PROBABILITIES |\n"
			"|___________________________|\n\n" );

	for( int j = 0; j < this->getLength(); j++ ){
		std::cout << "At position " << j+1 << std::endl;
		for( int k = 0; k < Global::modelOrder + 1; k++ ){
			for( int y = 0; y < pow( Alphabet::getSize(), k + 1 ); y++ ){
				std::cout << std::fixed << std::setprecision(4) << v_[k][y][j] << '\t';
			}
			std::cout << std::endl;
		}
	}
}

void Motif::write(){
	std::string opath = std::string( Global::outputDirectory )  + '/' + std::string( Global::posSequenceBasename ) + ".condsInit";
	std::ofstream ofile( opath.c_str() );
	for( int j = 0; j < length_; j++ ){
		for( int k = 0; k < Global::modelOrder + 1; k++ ){
			for( int y = 0; y < pow( Alphabet::getSize(), k+1 ); y++ ){
				ofile << std::fixed << std::setprecision(4) << v_[k][y][j] << '\t';
			}
			ofile << std::endl;
		}
		ofile << std::endl;
	}
}

