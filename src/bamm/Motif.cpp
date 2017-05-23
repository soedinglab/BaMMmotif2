#include <fstream>		// std::fstream

#include "Motif.h"

Motif::Motif( int length ){

	int k, y;

	for( k = 0; k < Global::Yk; k++ ){
		Y_.push_back( ipow( Alphabet::getSize(), k ) );
	}

	W_ = length;

	int K = Global::modelOrder;

	v_ = ( float*** )calloc( K+1, sizeof( float** ) );
	n_ = ( int*** )calloc( K+1, sizeof( int** ) );
	p_ =  ( float*** )calloc( K+1, sizeof( float** ) );
	for( k = 0; k < K+1; k++ ){
		v_[k] = ( float** )calloc( Y_[k+1], sizeof( float* ) );
		n_[k] = ( int** )calloc( Y_[k+1], sizeof( int* ) );
		p_[k] = ( float** )calloc( Y_[k+1], sizeof( float* ) );
		for( y = 0; y < Y_[k+1]; y++ ){
			v_[k][y] = ( float* )calloc( W_, sizeof( float ) );
			n_[k][y] = ( int* )calloc( W_+ K, sizeof( int ) ) + K;
			p_[k][y] = ( float* )calloc( W_, sizeof( float ) );
		}
	}

	s_ = ( float** )calloc( Y_[K+1]+1, sizeof( float* ) )+1;
	for( int y = -1; y < Y_[K+1]; y++ ){
		s_[y] = ( float* )calloc( W_, sizeof( float ) );
	}

	f_bg_ = Global::negSequenceSet->getBaseFrequencies();
}

Motif::Motif( const Motif& other ){ 		// copy constructor

	int k, y, j;

	for( k = 0; k < Global::Yk; k++ ){
		Y_.push_back( ipow( Alphabet::getSize(), k ) );
	}

	W_ = other.W_;

	int K = Global::modelOrder;

	v_ = ( float*** )malloc( ( K+1 ) * sizeof( float** ) );
	n_ = ( int*** )malloc( ( K+1 ) * sizeof( int** ) );
	p_ = ( float*** )malloc( ( K+1 )* sizeof( float** ) );
	for( k = 0; k < K+1; k++ ){
		v_[k] = ( float** )malloc( Y_[k+1] * sizeof( float* ) );
		n_[k] = ( int** )malloc( Y_[k+1] * sizeof( int* ) );
		p_[k] = ( float** )malloc( Y_[k+1] * sizeof( float* ) );
		for( y = 0; y < Y_[k+1]; y++ ){
			v_[k][y] = ( float* )malloc( W_ * sizeof( float ) );
			n_[k][y] = ( int* )malloc( ( W_ + K ) * sizeof( int ) ) + K;
			p_[k][y] = ( float* )malloc( W_ * sizeof( float ) );
			for( j = 0; j < W_; j++ ){
				v_[k][y][j] = other.v_[k][y][j];
				p_[k][y][j] = other.p_[k][y][j];
			}
			for( j = -K; j < W_; j++ ){
				n_[k][y][j] = other.n_[k][y][j];
			}
		}
	}

	s_ = ( float** )calloc( Y_[K+1]+1, sizeof( float* ) )+1;
	for( int y = -1; y < Y_[K+1]; y++ ){
		s_[y] = ( float* )calloc( W_, sizeof( float ) );
		for( j = 0; j < W_; j++ ){
			s_[y][j] = other.s_[y][j];
		}
	}

	f_bg_ = other.f_bg_;

	isInitialized_ = true;
}

Motif::~Motif(){

	for( int k = 0; k < Global::modelOrder+1; k++ ){
		for( int y = 0; y < Y_[k+1]; y++ ){
			free( v_[k][y] );
			free( n_[k][y] - Global::modelOrder );
			free( p_[k][y] );
		}
		free( v_[k] );
		free( n_[k] );
		free( p_[k] );
	}
	free( v_ );
	free( n_ );
	free( p_ );

	for( int y = 0; y < Y_[Global::modelOrder+1]; y++ ){
		free( s_[y] );
	}
	free( s_-1 );

}

// initialize v from IUPAC pattern (BaMM pattern)
void Motif::initFromBaMMPattern( std::string pattern ){

	// todo: under construction
	// calculate v from the kmers in the pattern
	// for k = 0:
	for( size_t j = 0; j < pattern.length(); j++ ){
		if( pattern[j] == 'A' ){
			v_[0][0][j] = 1.0f;
		} else if( pattern[j] == 'C' ){
			v_[0][1][j] = 1.0f;
		} else if( pattern[j] == 'G' ){
			v_[0][2][j] = 1.0f;
		} else if( pattern[j] == 'T' ){
			v_[0][3][j] = 1.0f;
		} else if( pattern[j] == 'S' ){
			v_[0][1][j] = 0.5f;
			v_[0][2][j] = 0.5f;
		} else if( pattern[j] == 'W' ){
			v_[0][0][j] = 0.5f;
			v_[0][3][j] = 0.5f;
		} else if( pattern[j] == 'R' ){
			v_[0][0][j] = 0.5f;
			v_[0][2][j] = 0.5f;
		} else if( pattern[j] == 'Y' ){
			v_[0][1][j] = 0.5f;
			v_[0][3][j] = 0.5f;
		} else if( pattern[j] == 'M' ){
			v_[0][0][j] = 0.5f;
			v_[0][1][j] = 0.5f;
		} else if( pattern[j] == 'K' ){
			v_[0][2][j] = 0.5f;
			v_[0][3][j] = 0.5f;
		}
	}
	// for k > 0:
	for( int k = 1; k < Global::modelOrder+1; k++ ){
		for( size_t j = k; j < pattern.length(); j++ ){

		}
	}

	// set isInitialized
	isInitialized_ = true;
}

// initialize v from binding sites file
void Motif::initFromBindingSites( char* filename ){

	std::ifstream file( filename );						// read file
	std::string bindingsite;							// read each binding site sequence from each line
	int bindingSiteWidth;								// length of binding site from each line
	int i, y, k, j, n;
	int minL = Global::posSequenceSet->getMinL();

	while( getline( file, bindingsite ).good() ){

		C_++;											// count the number of binding sites

		// add alphabets randomly at the beginning of each binding site
		for( i = 0; i < Global::addColumns.at(0); i++ )
			bindingsite.insert( bindingsite.begin(), Alphabet::getBase( static_cast<uint8_t>( rand() % Y_[1] + 1) ) );

		// add alphabets randomly at the end of each binding site
		for( i = 0; i < Global::addColumns.at(1); i++ )
			bindingsite.insert( bindingsite.end(), Alphabet::getBase( static_cast<uint8_t>( rand() % Y_[1] + 1) ) );

		bindingSiteWidth = static_cast<int>( bindingsite.length() );

		if( bindingSiteWidth != W_ ){					// all the binding sites should have the same length
			fprintf( stderr, "Error: Length of binding site on line %d differs.\n"
					"Binding sites should have the same length.\n", C_ );
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
				for( n = k; n >= 0; n-- ){				// calculate y based on (k+1)-mer bases
					y += Y_[n] * ( Alphabet::getCode( bindingsite[j-n] ) - 1 );
				}
				n_[k][y][j]++;
			}
		}
	}

	// calculate v and p from k-mer counts n
	calculateV();
	calculateP();

	// set isInitialized
	isInitialized_ = true;

	// optional: save initial model
	if( Global::saveInitialModel ){
		write( -1 );
	}
}

// initialize v from PWM file
void Motif::initFromPWM( float** PWM, int asize, int count ){

	int k, y, j, n, i;
	int K = Global::modelOrder;
	// for k = 0, obtain v from PWM:
	for( j = 0; j < W_; j++ ){
		float norm = 0.0f;
		for( y = 0; y < asize; y++ ){
			v_[0][y][j] = PWM[y][j];
			norm += PWM[y][j];
		}
		// normalize PWMs to sum up the weights to 1
		for( y = 0; y < asize; y++ ){
			v_[0][y][j] /= norm;
		}
	}

	// learn higher-order model based on the weights from the PWM
	// compute log odd scores s[y][j], log likelihoods of the highest order K
	std::vector<std::vector<float>> score;
	score.resize( asize );
	for( y = 0; y < asize; y++ ){
		score[y].resize( W_ );
	}
	for( y = 0; y < asize; y++ ){
		for( j = 0; j < W_; j++ ){
			score[y][j] = v_[0][y][j] / f_bg_[y];
		}
	}

	// sampling z from each sequence of the sequence set based on the weights:
	std::vector<Sequence*> posSet = Global::posSequenceSet->getSequences();
	int N = Global::posSequenceSet->getN();

	for( n = 0; n < N; n++ ){

		int LW1 = posSet[n]->getL() - W_ + 1;

		// motif position
		int z;

		// get the kmer array
		int* kmer = posSet[n]->getKmer();

		// calculate responsibilities over all LW1 positions on n'th sequence
		std::vector<float> r;
		r.resize( LW1 + 1 );

		std::vector<float> posteriors;
		float normFactor = 0.0f;
		// calculate positional prior:
		float pos0 = 1.0f - Global::q;
		float pos1 = Global::q / static_cast<float>( LW1 );

		// todo: could be parallelized by extracting 8 sequences at once
		// todo: should be written in a faster way
		for( i = 1; i <= LW1; i++ ){
			r[i] = 1.0f;
			for( j = 0; j < W_; j++ ){
				// extract monomers on the motif at position i over W of the n'th sequence
				y = ( kmer[i-1+j] >= 0 ) ? kmer[i-1+j] % asize : -1;
				if( y >= 0 )	r[i] *= score[y][j];
			}
			r[i] *= pos1;
			normFactor += r[i];
		}
		// for sequences that do not contain motif
		r[0] = pos0;
		normFactor += r[0];
		for( i = 0; i <= LW1; i++ ){
			r[i] /= normFactor;
			posteriors.push_back( r[i] );
		}
		// draw a new position z from discrete posterior distribution
		std::discrete_distribution<> posterior_dist( posteriors.begin(), posteriors.end() );

		// draw a sample z randomly
		z = posterior_dist( Global::rngx );

		// count kmers with sampled z
		if( z > 0 ){
			for( k = 0; k < K+1; k++ ){
				for( j = ( z <= K ) ? 1-z : -K; j < W_; j++ ){
					y = ( kmer[z-1+j] >= 0 ) ? kmer[z-1+j] % Y_[k+1] : -1;
					if( y >= 0 ){
						n_[k][y][j]++;
					}
				}
			}
		}
	}

	// calculate motif model from counts for higher order
	// for k > 0:
	for( k = 1; k < K+1; k++ ){
		for( y = 0; y < Y_[k+1]; y++ ){
			int y2 = y % Y_[k];									// cut off the first nucleotide in (k+1)-mer y
			int yk = y / Y_[1];									// cut off the last nucleotide in (k+1)-mer y
			for( j = 0; j < W_; j++ ){
				v_[k][y][j] = ( static_cast<float>( n_[k][y][j] ) + Global::modelAlpha.at(k) * v_[k-1][y2][j] )
							/ ( static_cast<float>( n_[k-1][yk][j-1] ) + Global::modelAlpha.at(k) );
			}
		}
	}

	// calculate probabilities p
	calculateP();

	// set isInitialized
	isInitialized_ = true;

}

// initialize v from Bayesian Markov model file and set isInitialized
void Motif::initFromBaMM( char* filename ){

	int k, y, j;
	std::ifstream file;
	file.open( filename, std::ifstream::in );
	std::string line;

	// loop over motif position j
	for( j = 0; j < W_; j++ ){

		// loop over order k
		for( k = 0; k < Global::modelOrder+1 ; k++ ){

			getline( file, line );

			std::stringstream number( line );

			float probability;

			y = 0;
			while( number >> probability ){
				// fill up for kmers y
				v_[k][y][j] = probability;
				//v_[k][y][j] = 1.0f / ( float )ipow( 4, k );
				y++;
			}
		}
		// read 'empty' line
		getline( file, line );
	}

	// calculate probabilities p
	calculateP();

	// set isInitialized
	isInitialized_ = true;

}

int Motif::getC(){
	return C_;
}

float** Motif::getS(){
	return s_;
}

float*** Motif::getP(){
	return p_;
}

int*** Motif::getN(){
	return n_;
}

void Motif::calculateV(){

	int y, j, k, y2, yk;

	// for k = 0, v_ = freqs:
	for( y = 0; y < Y_[1]; y++ ){
		for( j = 0; j < W_; j++ ){
			v_[0][y][j] = ( static_cast<float>( n_[0][y][j] ) + Global::modelAlpha.at(0) * f_bg_[y] )
						/ ( static_cast<float>( C_ ) + Global::modelAlpha.at(0) );
		}
	}

	// for k > 0:
	for( k = 1; k < Global::modelOrder+1; k++ ){
		for( y = 0; y < Y_[k+1]; y++ ){
			y2 = y % Y_[k];									// cut off the first nucleotide in (k+1)-mer y
			yk = y / Y_[1];									// cut off the last nucleotide in (k+1)-mer y
			for( j = 0; j < W_; j++ ){
				v_[k][y][j] = ( static_cast<float>( n_[k][y][j] ) + Global::modelAlpha.at(k) * v_[k-1][y2][j] )
							/ ( static_cast<float>( n_[k-1][yk][j-1] ) + Global::modelAlpha.at(k) );
			}
		}
	}
}

void Motif::calculateP(){
	// calculate probabilities, i.e. p(ACG) = p(G|AC) * p(AC)
	int j, y, k, y2, yk;

	// when k = 0:
	for( j = 0; j < W_; j++ ){
		for( y = 0; y < Y_[1]; y++ ){
			p_[0][y][j] = v_[0][y][j];
		}
	}
	// when k > 0:
	for( k = 1; k < Global::modelOrder+1; k++){
		for( y = 0; y < Y_[k+1]; y++ ){
			y2 = y % Y_[k];									// cut off the first nucleotide in (k+1)-mer
			yk = y / Y_[1];									// cut off the last nucleotide in (k+1)-mer
			for( j = 0; j < k; j++ ){
				p_[k][y][j] = p_[k-1][y2][j];				// i.e. p(ACG) = p(CG)
			}
			for( j = k; j < W_; j++ ){
				p_[k][y][j] =  v_[k][y][j] * p_[k-1][yk][j-1];
			}
		}
	}
}

void Motif::calculateS( float** Vbg ){

	int K = Global::modelOrder;
	int K_bg = ( Global::bgModelOrder < K ) ? Global::bgModelOrder : K ;

	for( int y = 0; y < Y_[K+1]; y++ ){
		int y_bg = y % Y_[K_bg+1];
		for( int j = 0; j < W_; j++ ){
			s_[y][j] = logf( v_[K][y][j] + 0.000001f ) - logf( Vbg[K_bg][y_bg] );
		}
	}

	for( int j = 0; j < W_; j++ ){
		s_[-1][j] = 0.0f;
	}

}

void Motif::calculateLinearS( float** Vbg ){

	int K = Global::modelOrder;
	int K_bg = ( Global::bgModelOrder < K ) ? Global::bgModelOrder : K ;

	for( int y = 0; y < Y_[K+1]; y++ ){
		int y_bg = y % Y_[K_bg+1];
		for( int j = 0; j < W_; j++ ){
			s_[y][j] =  v_[K][y][j] / Vbg[K_bg][y_bg];
		}
	}

	for( int j = 0; j < W_; j++ ){
		s_[-1][j] = 0.0f;
	}
}

void Motif::print(){

	fprintf( stderr," _______________________\n"
			"|                       |\n"
			"|  n for Initial Model  |\n"
			"|_______________________|\n\n" );
	for( int j = 0; j < W_; j++ ){
		for( int k = 0; k < Global::modelOrder+1; k++ ){
			for( int y = 0; y < Y_[k+1]; y++ )
				std::cout << std::scientific << n_[k][y][j] << '\t';
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
}

void Motif::write( int N ){						// write each motif with a number

	/**
	 * save motif learned by BaMM in two flat files:
	 * (1) posSequenceBasename.ihbcp: 		conditional probabilities after EM
	 * (2) posSequenceBasename.ihbp: 		probabilities of PWM after EM
	 */

	std::string opath = std::string( Global::outputDirectory )  + '/'
			+ Global::posSequenceBasename + "_motif_" + std::to_string( N+1 );

	// output conditional probabilities v[k][y][j] and probabilities prob[k][y][j]
	std::string opath_v = opath + ".ihbcp"; 	// inhomogeneous bamm conditional probabilities
	std::string opath_p = opath + ".ihbp";		// inhomogeneous bamm probabilities
	std::ofstream ofile_v( opath_v.c_str() );
	std::ofstream ofile_p( opath_p.c_str() );
	int j, k, y;

	std::vector<int> Y;
	for( k = 0; k < std::max( Global::modelOrder+2, Global::bgModelOrder+2 ); k++ ){
		Y.push_back( ipow( Alphabet::getSize(), k ) );
	}
	for( j = 0; j < W_; j++ ){
		for( k = 0; k < Global::modelOrder+1; k++ ){
			for( y = 0; y < Y[k+1]; y++ ){
				ofile_v << std::scientific << std::setprecision(3) << v_[k][y][j] << ' ';
				ofile_p << std::scientific << std::setprecision(3) << p_[k][y][j] << ' ';
			}
			ofile_v << std::endl;
			ofile_p << std::endl;
		}
		ofile_v << std::endl;
		ofile_p << std::endl;
	}
}
