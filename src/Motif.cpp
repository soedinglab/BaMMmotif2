/*
 * Motif.cpp
 *
 *  Created on: Apr 19, 2016
 *      Author: administrator
 */

#include "Motif.h"

Motif::Motif( int length ){

	int k, y;

	W_ = length;

	v_ = ( float*** )calloc( K_+1, sizeof( float** ) );
	n_ = ( int*** )calloc( K_+1, sizeof( int** ) );
	for( k = 0; k < K_+1; k++ ){
		v_[k] = ( float** )calloc( Global::powA[k+1], sizeof( float* ) );
		n_[k] = ( int** )calloc( Global::powA[k+1], sizeof( int* ) );
		for( y = 0; y < Global::powA[k+1]; y++ ){
			v_[k][y] = ( float* )calloc( W_, sizeof( float ) );
			n_[k][y] = ( int* )calloc( W_, sizeof( int ) );
		}
	}

	v_bg_ = ( float* )calloc( Global::powA[1], sizeof( float ) );

	BackgroundModel bg;
	for( y = 0; y < Global::powA[1]; y++ )
		v_bg_[y] = bg.getVbg()[0][y];
}

Motif::Motif( const Motif& other ){ 		// deep copy

	int k, y, j;

	W_ = other.W_;

	v_ = ( float*** )malloc( ( K_+1 )* sizeof( float** ) );
	n_ = ( int*** )malloc( ( K_+1 ) * sizeof( int** ) );
	for( k = 0; k < K_+1; k++ ){
		v_[k] = ( float** )malloc( Global::powA[k+1] * sizeof( float* ) );
		n_[k] = ( int** )malloc( Global::powA[k+1] * sizeof( int* ) );
		for( y = 0; y < Global::powA[k+1]; y++ ){
			v_[k][y] = ( float* )malloc( W_ * sizeof( float ) );
			n_[k][y] = ( int* )malloc( W_ * sizeof( int ) );
			for( j = 0; j < W_; j++ ){
				v_[k][y][j] = other.v_[k][y][j];
				n_[k][y][j] = other.n_[k][y][j];
			}
		}
	}

	isInitialized_ = true;

	v_bg_ = ( float* )calloc( Global::powA[1], sizeof( float ) );
	for( y = 0; y < Global::powA[1]; y++ )
		v_bg_[y] = other.v_bg_[y];
}

Motif::~Motif(){
	for( int k = 0; k < K_+1; k++ ){
		for( int y = 0; y < Global::powA[k+1]; y++ ){
			free( v_[k][y] );
			free( n_[k][y] );
		}
		free( v_[k] );
		free( n_[k] );
	}
	free( v_ );
	free( n_ );

	free( v_bg_ );
	std::cout << "Destructor for Motif class works fine. \n";
}

// initialize v from IUPAC pattern (BaMM pattern)
void Motif::initFromBaMMPattern( char* pattern ){
	// calculate k-mer counts n
	// calculate v from k-mer counts n using calculateV()
	// set isInitialized
}

// initialize v from binding sites file
void Motif::initFromBindingSites( char* filename ){

	std::ifstream file( filename );						// read file
	std::string bindingsite;							// read each binding site sequence from each line
	int bindingSiteWidth;								// length of binding site from each line
	int minL = Global::posSequenceSet->getMinL();
	int i, y, k, j, n;

	while( getline( file, bindingsite ).good() ){

		N_++;											// count the number of binding sites

		// add alphabets randomly at the beginning of each binding site
		for( i = 0; i < Global::addColumns.at(0); i++ )
			bindingsite.insert( bindingsite.begin(), Alphabet::getBase( rand() % Global::powA[1] + 1 ) );
		// add alphabets randomly at the end of each binding site
		for( i = 0; i < Global::addColumns.at(1); i++ )
			bindingsite.insert( bindingsite.end(), Alphabet::getBase( rand() % Global::powA[1] + 1 ) );

		bindingSiteWidth = bindingsite.length();

		if( bindingSiteWidth != W_ ){					// all the binding sites should have the same length
			fprintf( stderr, "Error: Length of binding site on line %d differs.\n"
					"Binding sites should have the same length.\n", N_ );
			exit( -1 );
		}
		if( bindingSiteWidth < K_+1 ){					// binding sites should be longer than the order of model
			fprintf( stderr, "Error: Length of binding site sequence "
					"is shorter than model order.\n" );
			exit( -1 );
		}
		if( bindingSiteWidth > minL ){					// binding sites should be shorter than the shortest posSeq
			fprintf( stderr, "Error: Length of binding site sequence "
					"exceeds the length of posSet sequence.\n" );
			exit( -1 );
		}

		// scan the binding sites and calculate k-mer counts n
		for( k = 0; k < K_+1; k++ ){					// k runs over all orders
			for( j = k; j < bindingSiteWidth; j++ ){	// j runs over all true motif positions
				y = 0;
				for( n = k; n >= 0; n-- )				// calculate y based on (k+1)-mer bases
					y += Global::powA[n] * ( Alphabet::getCode( bindingsite[j-n] ) - 1 );
				n_[k][y][j]++;
			}
		}
	}

	// calculate v from k-mer counts n
	calculateV();

	write();

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
	return W_;
}

void Motif::calculateV(){

	int y, j, k, y2, yk;
	// for k = 0, v_ = freqs:
	for( y = 0; y < Global::powA[1]; y++ )
		for( j = 0; j < W_; j++ )
			v_[0][y][j] = ( n_[0][y][j] + Global::modelAlpha.at(0) * v_bg_[y] )
						/ ( N_ + Global::modelAlpha.at(0) );

	// for k > 0:
	for( k = 1; k < K_+1; k++ ){
		for( y = 0; y < Global::powA[k+1]; y++ ){
			y2 = y % Global::powA[k];						// cut off the first nucleotide in (k+1)-mer y
			yk = y / Global::powA[1];						// cut off the last nucleotide in (k+1)-mer y
			for( j = 0; j < k; j++ )						// when j < k, i.e. p(A|CG) = p(A|C)
				v_[k][y][j] = v_[k-1][y2][j];
			for( j = k; j < W_; j++ )
				v_[k][y][j] = ( n_[k][y][j] + Global::modelAlpha.at(k) * v_[k-1][y2][j] )
							/ ( n_[k-1][yk][j-1] + Global::modelAlpha.at(k) );
		}
	}
}

float*** Motif::getV(){
	return v_;
}

// update v from fractional k-mer counts n and current alphas
void Motif::updateV( float*** n, float** alpha ){

	assert( isInitialized_ );

	int y, j, k, y2, yk;

	// sum up the n over (k+1)-mers at different position of motif
	float* sumN = ( float* )calloc( W_, sizeof( float ) );
	for( j = 0; j < W_; j++ )
		for( y = 0; y < Global::powA[1]; y++ )
			sumN[j] += n[0][y][j];

	// for k = 0, v_ = freqs:
	for( y = 0; y < Global::powA[1]; y++ ){
		for( j = 0; j < W_; j++ )
			v_[0][y][j] = ( n[0][y][j] + alpha[0][j] * v_bg_[y] )
						/ ( sumN[j] + alpha[0][j] );
	}

	free( sumN );

	// for k > 0:
	for( k = 1; k < K_+1; k++ ){
		for( y = 0; y < Global::powA[k+1]; y++ ){
			y2 = y % Global::powA[k];						// cut off the first nucleotide in (k+1)-mer
			yk = y / Global::powA[1];						// cut off the last nucleotide in (k+1)-mer
			for( j = 0; j < k; j++ )						// when j < k, i.e. p(A|CG) = p(A|C)
				v_[k][y][j] = v_[k-1][y2][j];
			for( j = k; j < W_; j++ )
				v_[k][y][j] = ( n[k][y][j] + alpha[k][j] * v_[k-1][y2][j] )
							/ ( n[k-1][yk][j-1] + alpha[k][j] );
		}
	}
}

void Motif::print(){
	printf( " _______________________\n"
			"|                       |\n"
			"|  v for Initial Model  |\n"
			"|_______________________|\n\n" );
	for( int j = 0; j < W_; j++ ){
		for( int k = 0; k < K_+1; k++ ){
			for( int y = 0; y < Global::powA[k+1]; y++ )
				std::cout << std::scientific << v_[k][y][j] << '\t';
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
}

void Motif::write(){

	/*
	 * save initial model in two flat files:
	 * (1) initialModelBasename.conds: 	conditional probabilities from initial model
	 * (2) initialModelBasename.counts:	counts of (k+1)-mers from initial model
	 */

	std::string opath = std::string( Global::outputDirectory )  + '/'
			+ std::string( Global::initialModelBasename );
	std::string opath_v = opath + ".conds";
	std::string opath_n = opath + ".counts";
	std::ofstream ofile_v( opath_v.c_str() );
	std::ofstream ofile_n( opath_n.c_str() );
	for( int j = 0; j < W_; j++ ){
		for( int k = 0; k < K_+1; k++ ){
			for( int y = 0; y < Global::powA[k+1]; y++ ){
				ofile_v << std::scientific << v_[k][y][j] << '\t';
				ofile_n << std::scientific << n_[k][y][j] << '\t';
			}
			ofile_v << std::endl;
			ofile_n << std::endl;
		}
		ofile_v << std::endl;
		ofile_n << std::endl;
	}
}
