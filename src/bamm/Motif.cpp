#include "Motif.h"

Motif::Motif( int length ){

	int k, y;
	int A = Alphabet::getSize();

	W_ = length;

	v_ = ( float*** )calloc( Global::modelOrder+1, sizeof( float** ) );
	n_ = ( int*** )calloc( Global::modelOrder+1, sizeof( int** ) );
	for( k = 0; k < Global::modelOrder+1; k++ ){
		v_[k] = ( float** )calloc( ipow( A, k+1 ), sizeof( float* ) );
		n_[k] = ( int** )calloc( ipow( A, k+1 ), sizeof( int* ) );
		for( y = 0; y < ipow( A, k+1 ); y++ ){
			v_[k][y] = ( float* )calloc( W_, sizeof( float ) );
			n_[k][y] = ( int* )calloc( W_, sizeof( int ) );
		}
	}

	v_bg_ = ( float* )calloc( A, sizeof( float ) );

	BackgroundModel bg( *Global::negSequenceSet,
						 Global::bgModelOrder,
						 Global::bgModelAlpha,
						 Global::interpolateBG );

	for( y = 0; y < A; y++ )
		v_bg_[y] = bg.getV()[0][y];

}

Motif::Motif( const Motif& other ){ 		// deep copy

	int k, y, j;
	int A = Alphabet::getSize();

	W_ = other.W_;

	v_ = ( float*** )malloc( ( Global::modelOrder+1 )* sizeof( float** ) );
	n_ = ( int*** )malloc( ( Global::modelOrder+1 ) * sizeof( int** ) );
	for( k = 0; k < Global::modelOrder+1; k++ ){
		v_[k] = ( float** )malloc( ipow( A, k+1 ) * sizeof( float* ) );
		n_[k] = ( int** )malloc( ipow( A, k+1 ) * sizeof( int* ) );
		for( y = 0; y < ipow( A, k+1 ); y++ ){
			v_[k][y] = ( float* )malloc( W_ * sizeof( float ) );
			n_[k][y] = ( int* )malloc( W_ * sizeof( int ) );
			for( j = 0; j < W_; j++ ){
				v_[k][y][j] = other.v_[k][y][j];
				n_[k][y][j] = other.n_[k][y][j];
			}
		}
	}

	isInitialized_ = true;

	v_bg_ = ( float* )calloc( A, sizeof( float ) );
	for( y = 0; y < A; y++ )
		v_bg_[y] = other.v_bg_[y];

}

Motif::~Motif(){

	int A = Alphabet::getSize();
	for( int k = 0; k < Global::modelOrder+1; k++ ){
		for( int y = 0; y < ipow( A, k+1 ); y++ ){
			free( v_[k][y] );
			free( n_[k][y] );
		}
		free( v_[k] );
		free( n_[k] );
	}
	free( v_ );
	free( n_ );

	free( v_bg_ );

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
	int i, y, k, j, n;
	int minL = Global::posSequenceSet->getMinL();
	int A = Alphabet::getSize();

	while( getline( file, bindingsite ).good() ){

		N_++;											// count the number of binding sites

		// add alphabets randomly at the beginning of each binding site
		for( i = 0; i < Global::addColumns.at(0); i++ )
			bindingsite.insert( bindingsite.begin(), Alphabet::getBase( static_cast<uint8_t>( rand() % A + 1) ) );

		// add alphabets randomly at the end of each binding site
		for( i = 0; i < Global::addColumns.at(1); i++ )
			bindingsite.insert( bindingsite.end(), Alphabet::getBase( static_cast<uint8_t>( rand() % A + 1) ) );

		bindingSiteWidth = static_cast<int>( bindingsite.length() );

		if( bindingSiteWidth != W_ ){					// all the binding sites should have the same length
			fprintf( stderr, "Error: Length of binding site on line %d differs.\n"
					"Binding sites should have the same length.\n", N_ );
			exit( -1 );
		}
		if( bindingSiteWidth < Global::modelOrder+1 ){	// binding sites should be longer than the order of model
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
		for( k = 0; k < Global::modelOrder+1; k++ ){	// k runs over all orders
			for( j = k; j < bindingSiteWidth; j++ ){	// j runs over all true motif positions
				y = 0;
				for( n = k; n >= 0; n-- )				// calculate y based on (k+1)-mer bases
					y += ipow( A, n ) * ( Alphabet::getCode( bindingsite[j-n] ) - 1 );
				n_[k][y][j]++;
			}
		}
	}

	// calculate v from k-mer counts n
	calculateV();

	// TODO: can be removed. By doing so, initial motif will not be output and latter will be replaced by refined model
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
	int A = Alphabet::getSize();

	// for k = 0, v_ = freqs:
	for( y = 0; y < A; y++ )
		for( j = 0; j < W_; j++ )
			v_[0][y][j] = ( static_cast<float>( n_[0][y][j] ) + Global::modelAlpha.at(0) * v_bg_[y] )
						/ ( static_cast<float>( N_ ) + Global::modelAlpha.at(0) );

	// for k > 0:
	for( k = 1; k < Global::modelOrder+1; k++ ){
		for( y = 0; y < ipow( A, k+1 ); y++ ){
			y2 = y % ipow( A, k );							// cut off the first nucleotide in (k+1)-mer y
			yk = y / A;										// cut off the last nucleotide in (k+1)-mer y
			for( j = 0; j < k; j++ )						// when j < k, i.e. p(A|CG) = p(A|C)
				v_[k][y][j] = v_[k-1][y2][j];
			for( j = k; j < W_; j++ )
				v_[k][y][j] = ( static_cast<float>( n_[k][y][j] ) + Global::modelAlpha.at(k) * v_[k-1][y2][j] )
							/ ( static_cast<float>( n_[k-1][yk][j-1] ) + Global::modelAlpha.at(k) );
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
	int A = Alphabet::getSize();

	// sum up the n over (k+1)-mers at different position of motif
	float* sumN = ( float* )calloc( W_, sizeof( float ) );
	for( j = 0; j < W_; j++ )
		for( y = 0; y < A; y++ )
			sumN[j] += n[0][y][j];

	// for k = 0, v_ = freqs:
	for( y = 0; y < A; y++ ){
		for( j = 0; j < W_; j++ )
			v_[0][y][j] = ( n[0][y][j] + alpha[0][j] * v_bg_[y] )
						/ ( sumN[j] + alpha[0][j] );
	}

	free( sumN );

	// for k > 0:
	for( k = 1; k < Global::modelOrder+1; k++ ){
		for( y = 0; y < ipow( A, k+1 ); y++ ){
			y2 = y % ipow( A, k );							// cut off the first nucleotide in (k+1)-mer
			yk = y / A;										// cut off the last nucleotide in (k+1)-mer
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
			"|  n for Initial Model  |\n"
			"|_______________________|\n\n" );
	for( int j = 0; j < W_; j++ ){
		for( int k = 0; k < Global::modelOrder+1; k++ ){
			for( int y = 0; y < ipow( Alphabet::getSize(), k+1 ); y++ )
				std::cout << std::scientific << n_[k][y][j] << '\t';
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
		for( int k = 0; k < Global::modelOrder+1; k++ ){
			for( int y = 0; y < ipow( Alphabet::getSize(), k+1 ); y++ ){
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
