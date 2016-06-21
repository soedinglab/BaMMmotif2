/*
 * Motif.cpp
 *
 *  Created on: Apr 19, 2016
 *      Author: administrator
 */

#include "Motif.h"

Motif::Motif( int length ){					// allocate memory for v_
	w_ = length;

	powA_ = new int[Global::modelOrder+2];
	for( int k = 0; k < Global::modelOrder + 2; k++ )
		powA_[k] = Global::ipow( Alphabet::getSize(), k );

	v_ = ( float*** )calloc( ( Global::modelOrder+1 ), sizeof( float** ) );
	for( int k = 0; k <= Global::modelOrder; k++ ){
		v_[k] = ( float** )calloc( powA_[k+1], sizeof( float* ) );
		for( int y = 0; y < powA_[k+1]; y++ )
			v_[k][y] = ( float* )calloc( w_, sizeof( float ) );
	}

	N_= 0;
}

Motif::Motif( const Motif& other ){ 		// deep copy
	w_ = other.w_;

	powA_ = new int[Global::modelOrder+2];
	for( int k = 0; k < Global::modelOrder + 2; k++ )
		powA_[k] = other.powA_[k];

	for( int k = 0; k < Global::modelOrder + 2; k++ )
		powA_[k] = Global::ipow( Alphabet::getSize(), k );
	if( other.v_ != NULL ){
		v_ = ( float*** )calloc( ( Global::modelOrder+1 ), sizeof( float** ) );
		for( int k = 0; k <= Global::modelOrder; k++ ){
			v_[k] = ( float** )calloc( powA_[k+1], sizeof( float* ) );
			for( int y = 0; y < powA_[k+1]; y++ ){
				v_[k][y] = ( float* )calloc( w_, sizeof( float ) );
				memcpy( v_[k][y], other.v_[k][y], w_ * sizeof( float ) );
			}
		}
		isInitialized_ = true;
	} else{
		v_ = NULL;
	}

	N_ = other.N_;
}

Motif::~Motif(){
	if( v_ ) free( v_ );
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
	int*** n;
	n = new int**[Global::modelOrder + 1];
	for( int k = 0; k < Global::modelOrder + 1; k++ ){
		n[k] = new int*[powA_[k+1]];
		for( int y = 0; y < powA_[k+1]; y++ ){
			n[k][y] = new int[w_];
			for( int j = 0; j < w_; j++ )
				n[k][y][j] = 0;
		}
	}

	std::ifstream file( filename );						// read file
	std::string bindingsite;							// read each binding site sequence from each line
	int w;												// length of binding site from each line

	while( getline( file, bindingsite ).good() ){

		N_++;											// count the number of binding sites

		// add alphabets randomly at the beginning of each binding site
		for( int i = 0; i < Global::addColumns.at(0); i++ )
			bindingsite.insert( bindingsite.begin(), Alphabet::getBase( rand() % powA_[1] + 1 ) );
		// add alphabets randomly at the end of each binding site
		for( int i = 0; i < Global::addColumns.at(1); i++ )
			bindingsite.insert( bindingsite.end(), Alphabet::getBase( rand() % powA_[1] + 1 ) );

		w = bindingsite.length();

		if( w != w_ ){									// all the binding sites should have the same length
			fprintf( stderr, "Error: Length of binding site sequence on line %d differs.\n"
					"Binding sites should have the same length.\n", N_ );
			exit( -1 );
		}
		if( w < Global::modelOrder + 1 ){				// binding sites should be longer than the order of model
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
		for( int k = 0; k < Global::modelOrder + 1; k++ ){				// k runs over all orders
			// j runs over all true motif positions
			for( int j = k; j < w; j++ ){
				int y = 0;
				for( int i = 0; i <= k; i++ )							// calculate y based on (k+1)-mer bases
					y += powA_[i] * ( Alphabet::getCode( bindingsite[j-k+i] ) - 1 );
				n[k][y][j]++;
			}
		}
	}

//	if( Global::verbose ){
//		printf( " ___________________________\n"
//				"|                           |\n"
//				"|  N_occ for Initial Model  |\n"
//				"|___________________________|\n\n" );
//		for( int k = 0; k < Global::modelOrder + 1; k++ ){
//			std::cout << "when k = " << k << ", n[k][y][j]"<< std::endl;
//			for( int y = 0; y < powA_[k+1]; y++ ){
//				for( int j = 0; j < w_; j++ )
//					std::cout << n[k][y][j] << '\t';
//				std::cout << std::endl;
//			}
//			std::cout << std::endl;
//		}
//	}

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

void Motif::calculateV( int*** n ){					// N is the number of motifs

	BackgroundModel bg;
	bg.init( std::vector<int> () );
	// for k = 0, v_ = freqs:
	for( int y = 0, k = 0; y < powA_[k+1]; y++ )
		for( int j = 0; j < w_; j++ )
			v_[k][y][j] = ( n[k][y][j] + Global::modelAlpha.at(k) * bg.getV()[y] ) / ( N_ + Global::modelAlpha.at(k) );

	// for k > 0:
	for( int k = 1; k < Global::modelOrder + 1; k++ ){
		for( int y = 0; y < powA_[k+1]; y++ ){
			int y2 = y / powA_[1];								// cut off the first nucleotide in (k+1)-mer y
			int yk = y % powA_[k];								// cut off the last nucleotide in (k+1)-mer y
			for( int j = 0; j < w_; j++ )
				v_[k][y][j] = ( n[k][y][j] + Global::modelAlpha.at(k) * v_[k-1][y2][j] ) / ( n[k-1][yk][j-1] + Global::modelAlpha.at(k) );
		}
	}

	// for normalization:
//	for( int k = 0; k < Global::modelOrder + 1; k++ )
//		for( int y = 0; y < powA_[k+1]; y++ )
//			for( int j = 0; j < w_; j++ )
//				v_[k][y][j] /= powA_[k] ;

//	if( Global::verbose ){
//		printf( " _______________________\n"
//				"|                       |\n"
//				"|  v for Initial Model  |\n"
//				"|_______________________|\n\n" );
//		for( int k = 0; k < Global::modelOrder + 1; k++ ){
//			std::cout << "when k = " << k << std::endl;
//			for( int y = 0; y < powA_[k+1]; y++ ){
//				for( int j = 0; j < w_; j++ )
//					std::cout << v_[k][y][j] << '\t';
//				std::cout << std::endl;
//			}
//			std::cout << std::endl;
//		}
//	}
}

float*** Motif::getV(){
	return v_;
}

void Motif::updateV( float*** n, float** alpha ){

	assert( isInitialized_ );

	// update v from fractional k-mer counts n and current alphas
	BackgroundModel bg;
	bg.init( std::vector<int> () );
	// for k = 0, v_ = freqs:
	for( int y = 0, k = 0; y < powA_[k+1]; y++ )
		for( int j = 0; j < w_; j++ )
			v_[k][y][j] = ( n[k][y][j] + alpha[k][j] * bg.getV()[y] ) / ( N_ + alpha[k][j] );

	// for k > 0:
	for( int k = 1; k < Global::modelOrder + 1; k++ ){
		for( int y = 0; y < powA_[k+1]; y++ ){
			int y2 = y / powA_[1];										// cut off the first nucleotide in (k+1)-mer
			int yk = y % powA_[k];										// cut off the last nucleotide in (k+1)-mer
			for( int j = 0; j < w_; j++ )
				v_[k][y][j] = ( n[k][y][j] + alpha[k][j] * v_[k-1][y2][j] ) / ( n[k-1][yk][j-1] + alpha[k][j] );
		}
	}
}

void Motif::print(){
	for( int k = 0; k < Global::modelOrder + 1; k++ ){
		std::cout << "when k = " << k << std::endl;
		for( int y = 0; y < powA_[k+1]; y++ ){
			for( int j = 0; j < w_; j++ )
				std::cout << v_[k][y][j] << '\t';
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	// only for testing:
//	for( int k = 0; k < Global::modelOrder + 1; k++ ){
//		std::cout << "when k = " << k << std::endl;
//		for( int j = 0; j < w_; j++ ){
//			float sum = 0.0f;
//			for( int y = 0; y < powA_[k+1]; y++ )
//				sum += v_[k][y][j];
//			std::cout << "at position " << j << ", sum = "<< sum << std::endl;
//		}
//	}
}

void Motif::write(){
	std::string opath = std::string( Global::outputDirectory )  + '/'
			+ std::string( Global::posSequenceBasename ) + ".condsInit";
	std::ofstream ofile( opath.c_str(), std::ios::app );	// append for MotifSet::write()
	for( int j = 0; j < w_; j++ ){
		for( int k = 0; k < Global::modelOrder + 1; k++ ){
			for( int y = 0; y < powA_[k+1]; y++ )
				ofile << std::scientific << v_[k][y][j] << '\t';
			ofile << std::endl;
		}
		ofile << std::endl;
	}
}
