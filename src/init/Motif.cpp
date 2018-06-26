#include <fstream>		// std::fstream
#include <random>       // std::mt19937, std::discrete_distribution
#include "Motif.h"

Motif::Motif( size_t length, size_t K, std::vector<float> alpha, float** v_bg, size_t k_bg, float glob_q ){

	W_ = length;
	K_ = K;
	C_ = 0;
    q_ = glob_q;

	for( size_t k = 0; k < K_+8; k++ ){
		Y_.push_back( ipow( Alphabet::getSize(), k ) );
	}

	if( v_bg != NULL ){
        k_bg_ = k_bg;
		v_bg_ = v_bg;
	} else {
        k_bg_ = 2;
		v_bg_ = ( float** )calloc( k_bg_ + 1, sizeof( float* ) );
		for( size_t k = 0; k <= k_bg_; k++ ){
            v_bg_[k] = ( float* )calloc( Y_[k+1], sizeof( float ) );
            for( size_t y = 0; y < Y_[k+1]; y++ ){
                v_bg_[k][y] = powf( 1.0f / Y_[1], (float)(k+1) );
            }
		}
	}

	v_ = ( float*** )calloc( K_+1, sizeof( float** ) );
	n_ = ( int*** )calloc( K_+1, sizeof( int** ) );
	p_ = ( float*** )calloc( K_+1, sizeof( float** ) );
	A_ = ( float** )calloc( K_+1, sizeof( float* ) );
	for( size_t k = 0; k < K_+1; k++ ){
		v_[k] = ( float** )calloc( Y_[k+1], sizeof( float* ) );
		n_[k] = ( int** )calloc( Y_[k+1], sizeof( int* ) );
		p_[k] = ( float** )calloc( Y_[k+1], sizeof( float* ) );
		for( size_t y = 0; y < Y_[k+1]; y++ ){
			v_[k][y] = ( float* )calloc( W_, sizeof( float ) );
			n_[k][y] = ( int* )calloc( W_, sizeof( int ) );
			p_[k][y] = ( float* )calloc( W_, sizeof( float ) );
		}
		A_[k] = ( float* )calloc( W_, sizeof( float ) );
		for( size_t j = 0; j < W_; j++ ){
			A_[k][j] = alpha[k];
		}
	}

	s_ = ( float** )calloc( Y_[K_+1], sizeof( float* ) );
	for( size_t y = 0; y < Y_[K_+1]; y++ ){
		s_[y] = ( float* )calloc( W_, sizeof( float ) );
	}

}

Motif::Motif( const Motif& other ){ 		// copy constructor

	W_ = other.W_;
	K_ = other.K_;
	A_ = other.A_;
	Y_ = other.Y_;
	C_ = other.C_;
    q_ = other.q_;

	v_ = ( float*** )malloc( ( K_+1 ) * sizeof( float** ) );
	n_ = ( int*** )malloc( ( K_+1 ) * sizeof( int** ) );
	p_ = ( float*** )malloc( ( K_+1 ) * sizeof( float** ) );
	A_ = ( float** )malloc( ( K_+1 ) * sizeof( float* ) );
	for( size_t k = 0; k < K_+1; k++ ){
		v_[k] = ( float** )malloc( Y_[k+1] * sizeof( float* ) );
		n_[k] = ( int** )malloc( Y_[k+1] * sizeof( int* ) );
		p_[k] = ( float** )malloc( Y_[k+1] * sizeof( float* ) );
		for( size_t y = 0; y < Y_[k+1]; y++ ){
			v_[k][y] = ( float* )malloc( W_ * sizeof( float ) );
			n_[k][y] = ( int* )malloc( W_ * sizeof( int ) );
			p_[k][y] = ( float* )malloc( W_ * sizeof( float ) );
			for( size_t j = 0; j < W_; j++ ){
				v_[k][y][j] = other.v_[k][y][j];
				p_[k][y][j] = other.p_[k][y][j];
				n_[k][y][j] = other.n_[k][y][j];
			}
		}
		A_[k] = ( float* )malloc( W_ * sizeof( float ) );
		for( size_t j = 0; j < W_; j++ ){
			A_[k][j] = other.A_[k][j];
		}
	}

	s_ = ( float** )malloc( Y_[K_+1] * sizeof( float* ) );
	for( size_t y = 0; y < Y_[K_+1]; y++ ){
		s_[y] = ( float* )malloc( W_ * sizeof( float ) );
		for( size_t j = 0; j < W_; j++ ){
			s_[y][j] = other.s_[y][j];
		}
	}

	k_bg_ = other.k_bg_;
    k_bg_ = other.k_bg_;
    for( size_t k = 0; k <= k_bg_; k++ ){
        v_bg_[k] = ( float* )calloc( Y_[k+1], sizeof( float ) );
        for( size_t y = 0; y < Y_[k+1]; y++ ){
            v_bg_[k][y] = other.v_bg_[k][y];
        }
    }
	isInitialized_ = true;
}

Motif::~Motif(){

	for( size_t k = 0; k < K_+1; k++ ){
		for( size_t y = 0; y < Y_[k+1]; y++ ){
			free( v_[k][y] );
			free( n_[k][y] );
			free( p_[k][y] );
		}
		free( v_[k] );
		free( n_[k] );
		free( p_[k] );
		free( A_[k]);
	}
	free( v_ );
	free( A_ );
	free( n_ );
	free( p_ );

	for( size_t y = 0; y < Y_[K_+1]; y++ ){
		free( s_[y] );
	}
	free( s_ );

/*  // does not work properly
    for( size_t k = 0; k <= k_bg_; k++ ){
        free( v_bg_[k] );
    }
    free( v_bg_ );*/
}

// initialize v from binding sites file
void Motif::initFromBindingSites( char* indir, size_t l_flank, size_t r_flank ){

	std::ifstream file( indir );			// read file
	std::string bindingsite;				// get binding site sequence from each line
	size_t bindingSiteWidth;				// length of binding site

	while( getline( file, bindingsite ).good() ){

		C_++;								// count the number of binding sites

		// add alphabets randomly at the beginning of each binding site
		for( size_t i = 0; i < l_flank; i++ )
			bindingsite.insert( bindingsite.begin(),
                                Alphabet::getBase( static_cast<uint8_t>( rand() )
                                                   % static_cast<uint8_t>( Y_[1] ) + 1 ) );

		// add alphabets randomly at the end of each binding site
		for( size_t i = 0; i < r_flank; i++ )
			bindingsite.insert( bindingsite.end(),
                                Alphabet::getBase( static_cast<uint8_t>( rand() )
                                                   % static_cast<uint8_t>( Y_[1] ) + 1 ) );

		bindingSiteWidth = bindingsite.length();

		if( bindingSiteWidth != W_ ){	// all binding sites have the same length
			fprintf( stderr, "Error: Length of binding site on line %d differs.\n"
					"Binding sites should have the same length.\n", ( int )C_ );
			exit( 1 );
		}
		if( bindingSiteWidth < K_+1 ){	// binding sites should be longer than the order of model
			fprintf( stderr, "Error: Length of binding site sequence "
					"is shorter than model order.\n" );
			exit( 1 );
		}

		// scan the binding sites and calculate k-mer counts n
		for( size_t k = 0; k < K_+1; k++ ){
			for( size_t j = k; j < bindingSiteWidth; j++ ){
				size_t y = 0;
				for( size_t a = 0; a < k+1; a++ ){
					// calculate y based on (k+1)-mer bases
					y += Y_[a] * ( Alphabet::getCode( bindingsite[j-a] ) - 1 );
				}
				n_[k][y][j]++;
			}
		}
	}

	// calculate v and p from k-mer counts n
	calculateV( n_ );
	calculateP();

	// set isInitialized
	isInitialized_ = true;

}

// initialize v from PWM file
void Motif::initFromPWM( float** PWM, size_t asize, SequenceSet* posSeqset, float q ){

    q_ = q;

	// set k-mer counts to zero
	for( size_t k = 0; k < K_+1; k++ ){
		for( size_t y = 0; y < Y_[k+1]; y++ ){
			for( size_t j = 0; j < W_; j++ ){
				n_[k][y][j] = 0;
			}
		}
	}

	// for k = 0, obtain v from PWM:
	for( size_t j = 0; j < W_; j++ ){
		float norm = 0.0f;
		for( size_t y = 0; y < asize; y++ ){
			if( PWM[y][j] <= 1.e-8 ){
                v_[0][y][j] = 1.e-8;
            } else {
                v_[0][y][j] = PWM[y][j];
            }
			norm += v_[0][y][j];
		}
		// normalize PWMs to sum the weights up to 1
		for( size_t y = 0; y < asize; y++ ){
			v_[0][y][j] /= norm;
		}
	}

	// learn higher-order model based on the weights from the PWM
	// compute log odd scores s[y][j], log likelihoods of the highest order K
	std::vector<std::vector<float>> score;
	score.resize( asize );
	for( size_t y = 0; y < asize; y++ ){
		score[y].resize( W_ );
	}
	for( size_t y = 0; y < asize; y++ ){
		for( size_t j = 0; j < W_; j++ ){
			score[y][j] = v_[0][y][j] / v_bg_[0][y];
		}
	}

	// sampling z from each sequence of the sequence set based on the weights:
	std::vector<Sequence*> posSet = posSeqset->getSequences();
	std::mt19937 rngx;
    size_t count = 0;
    std::vector<Sequence*>::iterator it = posSet.begin();
    while( it != posSet.end() ){
        if( (*it)->getL() < W_ ){
            //std::cout << "Warning: remove the short sequence: " << (*it)->getHeader() << std::endl;
            posSet.erase(it);
            count++;
        } else {
            it++;
        }
    }
    if( count > 0 ) std::cout << "Note: " << count
                              << " short sequences have been neglected for sampling PWM."
                              << std::endl;

#pragma omp parallel for

	for( size_t n = 0; n < posSet.size(); n++ ){

		size_t LW1 = posSet[n]->getL() - W_ + 1;

		// motif position
		size_t z;

		// get the kmer array
		size_t* kmer = posSet[n]->getKmer();

		// calculate responsibilities over all LW1 positions on n'th sequence
		std::vector<float> r;
		r.resize( LW1 + 1 );

		std::vector<float> posteriors;
		float normFactor = 0.0f;
		// calculate positional prior:
		float pos0 = 1.0f - q;
		float pos1 = q / static_cast<float>( LW1 );

		// todo: could be parallelized by extracting 8 sequences at once
		// todo: should be written in a faster way
		for( size_t i = 1; i <= LW1; i++ ){
			r[i] = 1.0f;
			for( size_t j = 0; j < W_; j++ ){
				// extract monomers from motif at position i
				// over W of the n'th sequence
				size_t y = kmer[i-1+j] % asize;
				r[i] *= score[y][j];
			}
			r[i] *= pos1;
			normFactor += r[i];
		}
		// for sequences that do not contain motif
		r[0] = pos0;
		normFactor += r[0];
		for( size_t i = 0; i <= LW1; i++ ){
			r[i] /= normFactor;
			posteriors.push_back( r[i] );
		}
		// draw a new position z from discrete posterior distribution
		std::discrete_distribution<size_t> posterior_dist( posteriors.begin(), posteriors.end() );

		// draw a sample z randomly
		z = posterior_dist( rngx );

		// count kmers with sampled z
		if( z > 0 ){
			for( size_t k = 0; k < K_+1; k++ ){
				for( size_t j = 0; j < W_; j++ ){
					size_t y = kmer[z-1+j] % Y_[k+1];
                    __sync_fetch_and_add(&(n_[k][y][j]), 1);
				}
			}
		}
	}

	// calculate motif model from counts for higher order
	// for k > 0:
	for( size_t k = 1; k < K_+1; k++ ){
		for( size_t y = 0; y < Y_[k+1]; y++ ){
			size_t y2 = y % Y_[k];			// cut off first nt in (k+1)-mer y
			size_t yk = y / Y_[1];			// cut off last nt in (k+1)-mer y
			for( size_t j = 0; j < k; j++ ){// when j < k, i.e. p(A|CG) = p(A|C)
				v_[k][y][j] = v_[k-1][y2][j];
			}
			for( size_t j = k; j < W_; j++ ){
				v_[k][y][j] = ( n_[k][y][j] + A_[k][j] * v_[k-1][y2][j] )
							/ ( n_[k-1][yk][j-1] + A_[k][j] );
			}
		}
	}
	// calculate probabilities p
	calculateP();

	// set isInitialized
	isInitialized_ = true;

}

// initialize v from Bayesian Markov model file and set isInitialized
void Motif::initFromBaMM( char* indir, size_t l_flank, size_t r_flank ){

	std::ifstream file( indir, std::ifstream::in );
	std::string line;
    if( file.is_open() ) {
/*    getline( file, line );
    std::stringstream eachline( line );
    size_t k = 0;     // count for
    while( line.length() != 0 ){
        std::cout << k << std::endl;
        getline( file, line );
        k++;
    }*/

        // loop over motif position j
        // set each v to 0.25f in the flanking region
        for (size_t j = 0; j < l_flank; j++) {
            for (size_t k = 0; k < K_ + 1; k++) {
                for (size_t y = 0; y < Y_[k + 1]; y++) {
                    v_[k][y][j] = 1.0f / static_cast<float>( Y_[1] );
                }
            }
        }

        // read in the v's from the bamm file for the core region
        for (size_t j = l_flank; j < W_ - r_flank; j++) {

            // loop over order k
            for (size_t k = 0; k < K_ + 1; k++) {

                getline(file, line);

                std::stringstream number(line);

                for (size_t y = 0; y < Y_[k + 1]; y++) {
                    number >> v_[k][y][j];
                }
            }
            // read 'empty' line
            getline(file, line);
        }

        // set each v to 0.25f in the flanking region
        for (size_t j = W_ - r_flank; j < W_; j++) {
            for (size_t k = 0; k < K_ + 1; k++) {
                for (size_t y = 0; y < Y_[k + 1]; y++) {
                    v_[k][y][j] = 1.0f / static_cast<float>( Y_[1] );
                }
            }
        }

        // calculate probabilities p
        calculateP();

        // set isInitialized
        isInitialized_ = true;

    } else {
        std::cerr << "Error: Input BaMM file cannot be opened!" << std::endl;
        exit( 1 );
    }
}

float** Motif::getS(){
	return s_;
}

void Motif::calculateV( int*** n ){
	// Note: This function is written for reading in binding site files.

    // for k = 0, v_ = freqs:
	for( size_t y = 0; y < Y_[1]; y++ ){
		for( size_t j = 0; j < W_; j++ ){
			v_[0][y][j] = ( n[0][y][j] + A_[0][j] * v_bg_[0][y] )
						/ ( static_cast<float>( C_ ) + A_[0][j] );
        }
	}

	// for k > 0:
	for( size_t k = 1; k < K_+1; k++ ){
		for( size_t y = 0; y < Y_[k+1]; y++ ){
			size_t y2 = y % Y_[k];				// cut off first nucleotide in (k+1)-mer y
			size_t yk = y / Y_[1];				// cut off last nucleotide in (k+1)-mer y
			for( size_t j = 0; j < k; j++ ){	// when j < k, i.e. p_j(A|CG) = p_j(A|C)
				v_[k][y][j] = v_[k-1][y2][j];
			}
			for( size_t j = k; j < W_; j++ ){
				v_[k][y][j] = ( n[k][y][j] + A_[k][j] * v_[k-1][y2][j] )
							/ ( n[k-1][yk][j-1] + A_[k][j] );
			}
		}
	}
}

void Motif::calculateP(){
	// calculate probabilities, i.e. p_j(ACG) = p_j(G|AC) * p_j-1(AC)

	// when k = 0:
	for( size_t j = 0; j < W_; j++ ){
		for( size_t y = 0; y < Y_[1]; y++ ){
			p_[0][y][j] = v_[0][y][j];
		}
	}

	// when k > 0:
	for( size_t k = 1; k < K_+1; k++ ){
		for( size_t y = 0; y < Y_[k+1]; y++ ){
			size_t yk = y / Y_[1];				// cut off last nucleotide in (k+1)-mer
			// when j < k:
			for( size_t j = 0; j < k; j++ ){
				p_[k][y][j] = 1;	            // i.e. p_j(ACG) = p_j(G|AC) x p_j-1(C|A) x p_j-2(A)
                for( size_t i = 0; i <= j; i++ ){
                    size_t yi = y / Y_[i];
                    p_[k][y][j] *= v_[k-i][yi][j-i];
                }
                for( size_t i = j+1; i <= k; i++ ){

                    if( (k - i) <= k_bg_ or k <= k_bg_ ){
                        size_t yi = y / Y_[i];
                        p_[k][y][j] *= v_bg_[k-i][yi];
                    } else {
                        size_t yi = y / Y_[1] % Y_[k_bg_+1];
                        p_[k][y][j] *= v_bg_[k_bg_][yi];
                    }
                }
			}
            // when j >= k:
			for( size_t j = k; j < W_; j++ ){
				p_[k][y][j] = v_[k][y][j] * p_[k-1][yk][j-1];
			}

		}
	}
}

void Motif::calculateLogS( float** Vbg, size_t K_bg ){

	for( size_t y = 0; y < Y_[K_+1]; y++ ){
		size_t y_bg = y % Y_[K_bg+1];
		for( size_t j = 0; j < W_; j++ ){
            // todo: randomize this value between
            //float rand = ((float) rand() / (RAND_MAX))/ 1e6f;
            float rand = 1e-5f;
			s_[y][j] = logf( v_[K_][y][j] + rand ) - logf( Vbg[K_bg][y_bg] );
		}
	}

}

void Motif::calculateLinearS( float** Vbg, size_t K_bg ){

	for( size_t y = 0; y < Y_[K_+1]; y++ ){
		size_t y_bg = y % Y_[K_bg+1];
		for( size_t j = 0; j < W_; j++ ){
			s_[y][j] = v_[K_][y][j] / Vbg[K_bg][y_bg];
		}
	}

}

void Motif::print(){

	fprintf( stderr," _______________________\n"
					"|                       |\n"
					"|  v for Initial Model  |\n"
					"|_______________________|\n\n" );
	for( size_t j = 0; j < W_; j++ ){
		for( size_t k = 0; k < K_+1; k++ ){
            float sum = 0.f;
			for( size_t y = 0; y < Y_[k+1]; y++ ) {
                std::cout << std::scientific << v_[k][y][j] << '\t';
                sum += v_[k][y][j];
            }
			std::cout << "\t sum = " << sum << std::endl;
		}
		std::cout << std::endl;
	}
}

void Motif::write( char* odir, std::string basename ){

	/**
	 * save motif learned by BaMM in two flat files:
	 * (1) posSequenceBasename.ihbcp: 		conditional probabilities after EM
	 * (2) posSequenceBasename.ihbp: 		probabilities of PWM after EM
	 */

	std::string opath = std::string( odir )  + '/' + basename;

	// output conditional probabilities v[k][y][j] and probabilities p[k][y][j]
	std::string opath_v = opath + ".ihbcp"; 	// ihbcp: inhomogeneous bamm
												// conditional probabilities
	std::string opath_p = opath + ".ihbp";		// ihbp: inhomogeneous bamm
												// probabilities
	std::ofstream ofile_v( opath_v.c_str() );
	std::ofstream ofile_p( opath_p.c_str() );

	for( size_t j = 0; j < W_; j++ ){
		for( size_t k = 0; k < K_+1; k++ ){
			for( size_t y = 0; y < Y_[k+1]; y++ ){
				ofile_v << std::scientific << std::setprecision(3)
						<< v_[k][y][j] << ' ';
				ofile_p << std::scientific << std::setprecision(3)
						<< p_[k][y][j] << ' ';
			}
			ofile_v << std::endl;
			ofile_p << std::endl;
		}
		ofile_v << std::endl;
		ofile_p << std::endl;
	}
}
