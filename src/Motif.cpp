/*
 * Motif.cpp
 *
 *  Created on: Apr 19, 2016
 *      Author: administrator
 */

#include "Motif.h"

Motif::Motif( int length ){					// allocate memory for v_
	w_ = length;

	v_ = ( float*** )calloc( ( K_+1 ), sizeof( float** ) );
	for( int k = 0; k <= K_; k++ ){
		v_[k] = ( float** )calloc( Global::powA[k+1], sizeof( float* ) );
		for( int y = 0; y < Global::powA[k+1]; y++ )
			v_[k][y] = ( float* )calloc( w_, sizeof( float ) );
	}

	vbg_ = ( float* )calloc( Global::powA[1], sizeof( float ) );
	BackgroundModel bg;
	bg.init( std::vector<int> () );
	for( int y = 0; y < Global::powA[1]; y++ )
		vbg_[y] = bg.getV()[0][y];
}

Motif::Motif( const Motif& other ){ 		// deep copy
	w_ = other.w_;

	if( other.v_ != NULL ){
		v_ = ( float*** )calloc( ( K_+1 ), sizeof( float** ) );
		for( int k = 0; k <= K_; k++ ){
			v_[k] = ( float** )calloc( Global::powA[k+1], sizeof( float* ) );
			for( int y = 0; y < Global::powA[k+1]; y++ ){
				v_[k][y] = ( float* )calloc( w_, sizeof( float ) );
				memcpy( v_[k][y], other.v_[k][y], w_ * sizeof( float ) );
			}
		}
		isInitialized_ = true;
	} else{
		v_ = NULL;
	}
	vbg_ = other.vbg_;
}

Motif::~Motif(){
	for( int k = 0; k <= K_; k++ ){
		for( int y = 0; y < Global::powA[k+1]; y++ ){
			free( v_[k][y] );
		}
		free( v_[k] );
	}
	free( v_ );
	free( vbg_ );
}

// initialize v from IUPAC pattern (BaMM pattern)
void Motif::initFromBaMMPattern( char* pattern ){
	// calculate k-mer counts n
	// calculate v from k-mer counts n using calculateV(n)
	// set isInitialized
}

// initialize v from binding sites file
void Motif::initFromBindingSites( char* filename ){

	// initialize n
	int*** n = ( int*** )calloc( ( K_+1 ), sizeof( int** ) );
	for( int k = 0; k < K_+1; k++ ){
		n[k] = ( int** )calloc( Global::powA[k+1], sizeof( int* ) );
		for( int y = 0; y < Global::powA[k+1]; y++ )
			n[k][y] = ( int* )calloc( w_, sizeof( int ) );
	}

	std::ifstream file( filename );						// read file
	std::string bindingsite;							// read each binding site sequence from each line
	int w;												// length of binding site from each line

	while( getline( file, bindingsite ).good() ){

		N_++;											// count the number of binding sites

		// add alphabets randomly at the beginning of each binding site
		for( int i = 0; i < Global::addColumns.at(0); i++ )
			bindingsite.insert( bindingsite.begin(), Alphabet::getBase( rand() % Global::powA[1] + 1 ) );
		// add alphabets randomly at the end of each binding site
		for( int i = 0; i < Global::addColumns.at(1); i++ )
			bindingsite.insert( bindingsite.end(), Alphabet::getBase( rand() % Global::powA[1] + 1 ) );

		w = bindingsite.length();

		if( w != w_ ){									// all the binding sites should have the same length
			fprintf( stderr, "Error: Length of binding site on line %d differs.\n"
					"Binding sites should have the same length.\n", N_ );
			exit( -1 );
		}
		if( w < K_+1 ){									// binding sites should be longer than the order of model
			fprintf( stderr, "Error: Length of binding site sequence "
					"is shorter than model order.\n" );
			exit( -1 );
		}
		if( (unsigned int)w > Global::posSequenceSet->getMinL() ){		// binding sites should be shorter than the shortest posSeq
			fprintf( stderr, "Error: Length of binding site sequence "
					"exceeds the length of posSet sequence.\n" );
			exit( -1 );
		}

		// scan the binding sites and calculate k-mer counts n
		int y;
		for( int k = 0; k < K_+1; k++ ){				// k runs over all orders
			for( int j = k; j < w; j++ ){				// j runs over all true motif positions
				y = 0;
				for( int n = k; n >= 0; n-- )			// calculate y based on (k+1)-mer bases
					y += Global::powA[n] * ( Alphabet::getCode( bindingsite[j-n] ) - 1 );
				n[k][y][j]++;
			}
		}
	}

	if( Global::verbose ){
		printf( " ___________________________\n"
				"|                           |\n"
				"|  N_occ for Initial Model  |\n"
				"|___________________________|\n\n" );
		for( int k = 0; k < K_+1; k++ ){
			std::cout << "when k = " << k << std::endl;
			for( int y = 0; y < Global::powA[k+1]; y++ ){
				for( int j = 0; j < w_; j++ )
					std::cout << n[k][y][j] << '\t';
				std::cout << std::endl;
			}
		}
	}

	// calculate v and s from k-mer counts n
	calculateV( n );

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

int Motif::getN(){
	return N_;
}

int Motif::getW(){
	return w_;
}

void Motif::calculateV( int*** n ){
	// for k = 0, v_ = freqs:
	for( int y = 0; y < Global::powA[1]; y++ )
		for( int j = 0; j < w_; j++ )
			v_[0][y][j] = ( n[0][y][j] + Global::modelAlpha.at(0) * vbg_[y] )
						/ ( N_ + Global::modelAlpha.at(0) );

	// for k > 0:
	for( int k = 1; k < K_+1; k++ ){
		for( int y = 0; y < Global::powA[k+1]; y++ ){
			int y2 = y % Global::powA[k];						// cut off the first nucleotide in (k+1)-mer y
			int yk = y / Global::powA[1];						// cut off the last nucleotide in (k+1)-mer y
			for( int j = 0; j < k; j++ )						// when j < k, i.e. p(A|CG) = p(A|C)
				v_[k][y][j] = v_[k-1][y2][j];
			for( int j = k; j < w_; j++ )
				v_[k][y][j] = ( n[k][y][j] + Global::modelAlpha.at(k) * v_[k-1][y2][j] )
							/ ( n[k-1][yk][j-1] + Global::modelAlpha.at(k) );
		}
	}

	if( Global::verbose ){
		printf( " _______________________\n"
				"|                       |\n"
				"|  v for Initial Model  |\n"
				"|_______________________|\n\n" );
		for( int j = 0; j < w_; j++ ){
			for( int k = 0; k < K_+1; k++ ){
				for( int y = 0; y < Global::powA[k+1]; y++ )
					std::cout << std::scientific << v_[k][y][j] << '\t';
				std::cout << std::endl;
			}
			std::cout << std::endl;
		}
	}
}

float*** Motif::getV(){
	return v_;
}

// update v from fractional k-mer counts n and current alphas
void Motif::updateV( float*** n, float** alpha ){
	assert( isInitialized_ );
	// for k = 0, v_ = freqs:
	for( int y = 0, k = 0; y < Global::powA[k+1]; y++ )
		for( int j = 0; j < w_; j++ ){
			v_[k][y][j] = ( n[k][y][j] + alpha[k][j] * vbg_[y] )
						/ ( N_ + alpha[k][j] );
		}

	// for k > 0:
	for( int k = 1; k < K_+1; k++ ){
		for( int y = 0; y < Global::powA[k+1]; y++ ){
			int y2 = y % Global::powA[k];						// cut off the first nucleotide in (k+1)-mer
			int yk = y / Global::powA[1];						// cut off the last nucleotide in (k+1)-mer
			for( int j = 0; j < k; j++ )				// when j < k, i.e. p(A|CG) = p(A|C)
				v_[k][y][j] = v_[k-1][y2][j];
			for( int j = k; j < w_; j++ )
				v_[k][y][j] = ( n[k][y][j] + alpha[k][j] * v_[k-1][y2][j] )
							/ ( n[k-1][yk][j-1] + alpha[k][j] );
		}
	}
}

void Motif::print(){
	for( int j = 0; j < w_; j++ ){
		for( int k = 0; k < K_+1; k++ ){
			for( int y = 0; y < Global::powA[k+1]; y++ )
				std::cout << std::scientific << v_[k][y][j] << '\t';
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
}

void Motif::write(){
	std::string opath = std::string( Global::outputDirectory )  + '/'
			+ std::string( Global::posSequenceBasename ) + ".conds";
	std::ofstream ofile( opath.c_str(), std::ios::app );	// append for MotifSet::write()
	for( int j = 0; j < w_; j++ ){
		for( int k = 0; k < K_+1; k++ ){
			for( int y = 0; y < Global::powA[k+1]; y++ )
				ofile << std::scientific << v_[k][y][j] << '\t';
			ofile << std::endl;
		}
		ofile << std::endl;
	}
}
