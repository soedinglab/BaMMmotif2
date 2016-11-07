#include "CGS.h"

CGS::CGS( Motif* motif, BackgroundModel* bg, std::vector<int> folds ){

	motif_ = motif;
	bg_ = bg;

	for( int k = 0; k < std::max( Global::modelOrder+2,  Global::bgModelOrder+2 ); k++ ){
		Y_.push_back( ipow( Alphabet::getSize(), k ) );
	}

	int W = motif_->getW();
	int LW1 = Global::posSequenceSet->getMaxL() - W + 1;
	int N = Global::posSequenceSet->getN();

	// allocate memory for counts nz_[k][y][j]
	nz_ = ( int*** )calloc( Global::modelOrder+1, sizeof( int** ) );
	for( int k = 0; k < Global::modelOrder+1; k++ ){
		nz_[k] = ( int** )calloc( Y_[k+1], sizeof( int* ) );
		for( int y = 0; y < Y_[k+1]; y++ ){
			nz_[k][y] = ( int* )calloc( W, sizeof( int ) );
			for( int j = 0; j < W; j++ ){
				nz_[k][y][j] = motif_->getN()[k][y][j];
			}
		}
	}

	// allocate memory and initialize responsibility r_[n][i]
	r_ = ( float** )malloc( N * sizeof( float* ) );
	for( int n = 0; n < N; n++ ){
		r_[n] = ( float* )malloc( LW1 * sizeof( float ) );
		for( int i = 0; i < LW1; i++ ){
			r_[n][i] = 1.0f;
		}
	}

	// allocate memory and initialize alpha_[k][j]
	alpha_ = ( float** )malloc( ( Global::modelOrder+1 ) * sizeof( float* ) );
	for( int k = 0; k < Global::modelOrder+1; k++ ){
		alpha_[k] = ( float* )malloc( W * sizeof( float ) );
		for( int j = 0; j < W; j++ ){
			alpha_[k][j] = Global::modelAlpha[k];
		}
	}

	// allocate memory and initialize positional prior pos_[i]
	// TODO: this only fits to the case when sequences have the same length
	pos_ = ( float* )malloc( ( LW1+1 ) * sizeof( float ) );
	pos_[0] = 1 - q_;
	for( int i = 1; i <= LW1; i++ ){
		pos_[i] = q_ / ( float )LW1;
	}

}

CGS::~CGS(){

	for( int k = 0; k < Global::modelOrder+1; k++ ){
		for( int y = 0; y < Y_[k+1]; y++ ){
			free( nz_[k][y] );
		}
		free( nz_[k] );
		free( alpha_[k] );
	}
	free( nz_ );
	free( alpha_ );

	for( int n = 0; n < Global::posSequenceSet->getN(); n++ ){
		free( r_[n] );
	}
	free( r_ );

	free( pos_ );

}

void CGS::GibbsSampling(){

	printf( " ___________________________\n"
			"|                           |\n"
			"|  Collapsed Gibbs sampler  |\n"
			"|___________________________|\n\n" );

	clock_t t0 = clock();
	bool iterate = true;								// flag for iterating before convergence
	int N = Global::posSequenceSet->getN();
	std::vector<Sequence*> posSeqs = Global::posSequenceSet->getSequences();
	int W = motif_->getW();
	int K = Global::modelOrder;
	int K_bg = Global::bgModelOrder;
	int LW1 = Global::posSequenceSet->getMaxL()-W+1;	// TODO: This means Gibbs sampling only works for sequences with the same length!
//	float pos_0 = 1 - q_;								// positional prior when no motif is present on the sequence
	float pos_i = q_ / ( float )LW1; 					// positional prior when motif(s) is present on the sequence
	float v_diff;										// model parameter difference before and after iteration
	float** v_before;									// model parameter with highest order before iteration
	// allocate memory for parameters v[y][j] with the highest order
	v_before = ( float** )calloc( Y_[K+1], sizeof( float* ) );
	for( int y = 0; y < Y_[K+1]; y++ ){
		v_before[y] = ( float* )calloc( W, sizeof( float ) );
	}

	// iterate over
	while( iterate && ( CGSIterations_ < Global::maxCGSIterations ) ){

		CGSIterations_++;

		// get model parameter with highest order before sampling
		for( int y = 0; y < Y_[K+1]; y++ ){
			for( int j = 0; j < W; j++ ){
				v_before[y][j] = motif_->getV()[K][y][j];
			}
		}

		// choose a random position for z
		for( int z = 0; z < LW1; z++ ){

			// pick one sequence and update r
			for( int n = 0; n < N; n++ ){

				// initialize nz[K][y][j] with motif->getN()
				for( int y = 0; y < Y_[K+1]; y++ ){
					for( int j = 0; j < W; j++ ){
						nz_[K][y][j] = motif_->getN()[K][y][j];
					}
				}

				// calculate nz by excluding K-mer counts over W at position i of n'th sequence
				for( int j = 0; j < W; j++ ){
					// extract kmers on the motif at position i over W of the n'th sequence
					int y = posSeqs[n]->extractKmer( z+j, std::min( z+j, K ) );
					if( y >= 0 ){
						nz_[K][y][j]--;
					}
				}

				// calculate nz for lower order k
				for( int k = K; k > 0; k-- ){				// k runs over all lower orders
					for( int y = 0; y < Y_[k+1]; y++ ){
						int y2 = y % Y_[k];					// cut off the first nucleotide in (k+1)-mer
						for( int j = 0; j < W; j++ ){
							nz_[k-1][y2][j] += nz_[k][y][j];
						}
					}
				}

				// updated model parameters v excluding the n'th sequence
				motif_->updateV( nz_, alpha_ );

				// sampling equation: calculate responsibilities over all LW1 positions on n'th sequence
				for( int j = 0; j < W; j++ ){
					// extract kmers on the motif at position i over W of the n'th sequence
					int y = posSeqs[n]->extractKmer( z+j, std::min( z+j, K ) );
					r_[n][z] *= ( motif_->getV()[K][y][j] / bg_->getV()[std::min( K, K_bg )][y] );
				}
				r_[n][z] *= pos_i;
			}

			// check model parameter difference for convergence
			v_diff = 0.0f;
			for( int y = 0; y < Y_[K+1]; y++ ){
				for( int j = 0; j < W; j++ ){
					v_diff += fabsf( motif_->getV()[K][y][j] - v_before[y][j] );
				}
			}
			if( Global::verbose ){
				std::cout << CGSIterations_ << " iteration:	";
				std::cout << "at position " << z;
				std::cout << ", para_diff = " << v_diff << std::endl;
			}

			if( v_diff < Global::epsilon )		iterate = false;
		}


		// * optional: optimize hyper-parameter alpha
		if( !Global::noAlphaSampling )	alphaSampling();

		// * optional: sampling of hyper-parameter q
		if( !Global::noQSampling )		qSampling();

		// check model parameter difference for convergence
		v_diff = 0.0f;
		for( int y = 0; y < Y_[K+1]; y++ ){
			for( int j = 0; j < W; j++ ){
				v_diff += fabsf( motif_->getV()[K][y][j] - v_before[y][j] );
			}
		}

		if( Global::verbose ){
			std::cout << CGSIterations_ << " iteration:	";
			std::cout << "para_diff = " << v_diff << std::endl;
		}

		if( v_diff < Global::epsilon )		iterate = false;
	}
	fprintf( stdout, "\n--- Runtime for Collapsed Gibbs sampling: %.4f seconds ---\n", ( ( float )( clock() - t0 ) ) / CLOCKS_PER_SEC );
}


void CGS::alphaSampling(){

}

void CGS::qSampling(){

}

void CGS::print(){

}

void CGS::write(){

	/**
	 * save CGS parameters in three flat files:
	 * (1) posSequenceBasename.CGScounts:		refined fractional counts of (k+1)-mers
	 * (2) posSequenceBasename.CGSposterior:	responsibilities, posterior distributions
	 * (3) posSequenceBasename.CGSalpha:		hyper-parameter alphas
	 */

	int W = motif_->getW();
	std::vector<Sequence*> posSeqs = Global::posSequenceSet->getSequences();

	std::string opath = std::string( Global::outputDirectory ) + '/'
						+ std::string( Global::posSequenceBasename );

	// output (k+1)-mer counts nz[k][y][j]
	std::string opath_n = opath + ".CGScounts";
	std::ofstream ofile_n( opath_n.c_str() );
	for( int j = 0; j < W; j++ ){
		for( int k = 0; k < Global::modelOrder+1; k++ ){
			for( int y = 0; y < Y_[k+1]; y++ ){
				ofile_n << std::scientific << nz_[k][y][j] << ' ';
			}
			ofile_n << std::endl;
		}
		ofile_n << std::endl;
	}

	// output responsibilities r[n][i]
	std::string opath_r = opath + ".CGSposterior";
	std::ofstream ofile_r( opath_r.c_str() );
	for( int n = 0; n < Global::posSequenceSet->getN(); n++ ){
		for( int i = 0; i < posSeqs[n]->getL()-W+1; i++ ){
			ofile_r << std::scientific << r_[n][i] << ' ';
		}
		ofile_r << std::endl;
	}

	// output parameter alphas alpha[k][j]
	std::string opath_alpha = opath + ".CGSalpha";
	std::ofstream ofile_alpha( opath_alpha.c_str() );
	for( int k = 0; k < Global::modelOrder+1; k++ ){
		ofile_alpha << "k = " << k << std::endl;
		for( int j = 0; j < W; j++ ){
			ofile_alpha << std::setprecision( 3 ) << alpha_[k][j] << ' ';
		}
		ofile_alpha << std::endl;
	}

}

