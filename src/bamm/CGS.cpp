#include "CGS.h"
#include "EM.h"

#include <random>			// std::discrete_distribution

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
	n_z_ = ( int*** )calloc( Global::modelOrder+1, sizeof( int** ) );
	for( int k = 0; k < Global::modelOrder+1; k++ ){
		n_z_[k] = ( int** )calloc( Y_[k+1], sizeof( int* ) );
		for( int y = 0; y < Y_[k+1]; y++ ){
			n_z_[k][y] = ( int* )calloc( W, sizeof( int ) );
		}
	}

	// allocate memory and initialize responsibility r_[n][i]
	r_ = ( float** )calloc( N, sizeof( float* ) );
	for( int n = 0; n < N; n++ ){
		r_[n] = ( float* )calloc( LW1+1, sizeof( float ) );
	}

	// allocate memory for s_[y][j] and initialize it
	s_ = ( float** )calloc( Y_[Global::modelOrder+1], sizeof( float* ) );
	for( int y = 0; y < Y_[Global::modelOrder+1]; y++ ){
		s_[y] = ( float* )calloc( W, sizeof( float ) );
	}

	// allocate memory and initialize alpha_[k][j]
	alpha_ = ( float** )malloc( ( Global::modelOrder+1 ) * sizeof( float* ) );
	for( int k = 0; k < Global::modelOrder+1; k++ ){
		alpha_[k] = ( float* )malloc( W * sizeof( float ) );
		for( int j = 0; j < W; j++ ){
			alpha_[k][j] = Global::modelAlpha[k];
		}
	}

	// allocate memory for positional prior pos_[i]
	pos_ = ( float* )calloc( LW1+1, sizeof( float ) );

	// allocate memory for motif position z_[n]
	z_ = ( int* )calloc( N, sizeof( int ) );

}

CGS::~CGS(){

	for( int k = 0; k < Global::modelOrder+1; k++ ){
		for( int y = 0; y < Y_[k+1]; y++ ){
			free( n_z_[k][y] );
		}
		free( n_z_[k] );
		free( alpha_[k] );
	}
	free( n_z_ );
	free( alpha_ );

	for( int n = 0; n < Global::posSequenceSet->getN(); n++ ){
		free( r_[n] );
	}
	free( r_ );
	for( int y = 0; y < Y_[Global::modelOrder+1]; y++ ){
		free( s_[y] );
	}
	free( s_ );

	free( pos_ );
	free( z_ );

}

void CGS::GibbsSampling(){

	printf( " ___________________________\n"
			"|                           |\n"
			"|  Collapsed Gibbs sampler  |\n"
			"|___________________________|\n\n" );

	clock_t t0 = clock();
	bool iterate = true;								// flag for iterating before convergence

	EM em( motif_, bg_ );
	em.learnMotif();
	for( int n = 0; n < Global::posSequenceSet->getN(); n++ ){
//		z_[n] = 30;
		z_[n] = em.getZ()[n]+1;

	}

	// iterate over
	while( iterate && CGSIterations_ < Global::maxCGSIterations ){

		CGSIterations_++;

		if( Global::verbose /* && CGSIterations_ % 10 == 0*/ ){
			std::cout << std::endl << CGSIterations_ << " iteration:\t";
		}

		// sampling z and q
		sampling_z_q( CGSIterations_ );

		// only for writing out model after each iteration:
//		motif_->calculateP();
//		motif_->write( CGSIterations_ );

		// * optional: optimize hyper-parameter alpha
		if( !Global::noAlphaUpdating )	updateAlphas();

	}

	// calculate probabilities
	motif_->calculateP();

	fprintf( stdout, "\n--- Runtime for Collapsed Gibbs sampling: %.4f seconds ---\n", ( ( float )( clock() - t0 ) ) / CLOCKS_PER_SEC );
}

void CGS::sampling_z_q( int iteration ){

	int N = Global::posSequenceSet->getN();
	std::vector<Sequence*> posSeqs = Global::posSequenceSet->getSequences();
	int W = motif_->getW();
	int K = Global::modelOrder;
	int K_bg = Global::bgModelOrder;
	int n, k, y, y_prev, y2, y_bg, i, j, LW1;
	int N_0 = 0;								// counts of sequences which do not contain motifs.


/*	// reset n_z_[k][y][j]
	for( k = 0; k < K+1; k++ ){
		for( y = 0; y < Y_[k+1]; y++ ){
			for( j = 0; j < W; j++ ){
				n_z_[k][y][j] = 0;
			}
		}
	}
	for( i = 0; i < N; i++ ){
		for( j = 0; j < W; j++ ){
			y = posSeqs[i]->extractKmer( z_[i]+j, std::min( z_[i]+j, K ) );
			n_z_[K][y][j]++;
		}
	}
	// calculate nz for lower order k
	for( k = K; k > 0; k-- ){				// k runs over all lower orders
		for( y = 0; y < Y_[k+1]; y++ ){
			y2 = y % Y_[k];					// cut off the first nucleotide in (k+1)-mer
			for( j = 0; j < W; j++ ){
				n_z_[k-1][y2][j] += n_z_[k][y][j];
			}
		}
	}
	// print kmer counts out
	std::cout << std::endl;
	for( j = 0; j < W; j++ ){
	for( k = 0; k < K+1; k++ ){
		for( y = 0; y < Y_[k+1]; y++ ){

				std::cout << n_z_[k][y][j] << '\t';
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;*/

	// sampling z:
	// loop over all sequences and drop one sequence each time and update r
	for( n = 0; n < N; n++ ){
		LW1 = posSeqs[n]->getL() - W + 1;

		// calculate positional prior:
		pos_[0] =  ( float )( 1.0f - q_ );
		for( i = 1; i <= LW1; i++ ){
			pos_[i] = ( float )q_ / ( float )LW1;
		}
		// count K-mers at position z[i]+j except the n'th sequence
		if( n == 0 ){							// for the first sequence
			// reset n_z_[k][y][j] when K = 0
			for( y = 0; y < Y_[K+1]; y++ ){
				for( j = 0; j < W; j++ ){
					n_z_[K][y][j] = 0;
				}
			}
			for( i = 1; i < N; i++ ){
				for( j = 0; j < W; j++ ){
					y = posSeqs[i]->extractKmer( z_[i]-1+j, std::min( z_[i]-1+j, K ) );
					n_z_[K][y][j]++;
				}
			}
		} else {								// for the rest sequences
			for( j = 0; j < W; j++ ){
				if( z_[n-1] != 0 ){
					// add the kmer counts for the previous sequence with updated z
					y_prev = posSeqs[n-1]->extractKmer( z_[n-1]-1+j, std::min( z_[n-1]-1+j, K ) );
					n_z_[K][y_prev][j]++;
				}
				if( z_[n] != 0 ){
					// remove the kmer counts for the current sequence with old z
					y = posSeqs[n]->extractKmer( z_[n]-1+j, std::min( z_[n]-1+j, K ) );
					n_z_[K][y][j]--;
				}
			}
		}

		// reset n_z_[k][y][j] when k < K
		for( k = 0; k < K; k++ ){
			for( y = 0; y < Y_[k+1]; y++ ){
				for( j = 0; j < W; j++ ){
					n_z_[k][y][j] = 0;
				}
			}
		}
		// calculate nz for lower order k
		for( k = K; k > 0; k-- ){				// k runs over all lower orders
			for( y = 0; y < Y_[k+1]; y++ ){
				y2 = y % Y_[k];					// cut off the first nucleotide in (k+1)-mer
				for( j = 0; j < W; j++ ){
					n_z_[k-1][y2][j] += n_z_[k][y][j];
				}
			}
		}

		// updated model parameters v excluding the n'th sequence
		motif_->updateV( n_z_, alpha_ );

		// compute log odd scores s[y][j], log likelihoods of the highest order K
		for( y = 0; y < Y_[K+1]; y++ ){
			for( j = 0; j < W; j++ ){
				y_bg = y % Y_[K_bg+1];
				s_[y][j] = motif_->getV()[K][y][j] / bg_->getV()[std::min( K, K_bg )][y_bg];
			}
		}

		// sampling equation: calculate responsibilities over all LW1 positions on n'th sequence
		std::vector<float> posterior_array;
		float normFactor = 0.0f;
		for( i = 1; i <= LW1; i++ ){
			r_[n][i] = 1.0f;
			for( j = 0; j < W; j++ ){
				// extract k-mers on the motif at position i over W of the n'th sequence
				y = posSeqs[n]->extractKmer( i-1+j, std::min( i-1+j, K ) );
				r_[n][i] *= s_[y][j];
			}
			r_[n][i] *= pos_[i];
			normFactor += r_[n][i];
		}

		// for sequences that do not contain motif
		r_[n][0] = pos_[0];
		normFactor += r_[n][0];

		for( i = 0; i <= LW1; i++ ){
			r_[n][i] /= normFactor;
			posterior_array.push_back( r_[n][i] );
		}

		// draw a new position z from discrete posterior distribution
		std::discrete_distribution<> posterior_dist( posterior_array.begin(), posterior_array.end() );
		std::default_random_engine rand;		// pick a random number
		z_[n] = posterior_dist( rand );			// draw a sample z

		if( z_[n] == 0 ) N_0++;
	}

	// sampling q:
	// draw two random numbers Q and P from Gamma distribution
	std::gamma_distribution<> P_Gamma_dist( N_0 + 1, 1 );
	std::gamma_distribution<> Q_Gamma_dist( N - N_0 + 1, 1 );
	std::default_random_engine rand1;		// pick a random number
	std::default_random_engine rand2;		// pick another random number
	double P = P_Gamma_dist( rand1 );		// draw a sample for P
	double Q = Q_Gamma_dist( rand2 );		// draw a sample for Q

	q_ = Q / ( Q + P );						// calculate q_

	if( Global::verbose ){
		// checking z values from the first 20 sequences
		for( n = 0; n < 20; n++ ) std::cout << z_[n] << '\t';
		std::cout << N_0 << " sequences do not have motif. q = " << q_  << "\t";
		// only for testing:
/*		std::cout << std::endl << "nz[2][TGA][j] = ";
		for( j = 0; j < W; j++ ){
			std::cout << n_z_[K][56][j] << '\t';
		}*/

/*		std::cout << "\n After:";
		for( n = 0; n < N; n++ ){
			std::cout << "\n z[" << n << "] = " << z_[n] <<'\t';
			for( i = 0; i < LW1; i++ ){
				std::cout << "r["<<n<<"]["<<i <<"] = " << r_[n][i] <<'\t';
			}
			std::cout << '\n';
		}

		for( j = 0; j < W; j++ ){
			for( y = 0; y < Y_[K+1]; y++ ){

				std::cout << n_z_[K][y][j] << '\t';
			}
			std::cout << '\n';
		}
		std::cout << '\n';*/

	}
}

void CGS::updateAlphas(){
	// update alphas due to the learning rate eta and gradient of the log posterior of alphas

//	int K = Global::modelOrder;
//	int W = motif_->getW();

	// calcGrad_logPostAlphas();

}

float CGS::calcGrad_logPostAlphas( float alpha, int k, int j ){
	// calculate gradient of the log posterior of alphas
	float gradient_logPostAlphas;
	float*** v = motif_->getV();

	// the first term
	gradient_logPostAlphas = ( -2.0f ) / alpha;
	// the second term
	gradient_logPostAlphas += Global::modelBeta * powf( Global::modelGamma, ( float )k ) / powf( alpha, 2.0f );
	// the third term
	gradient_logPostAlphas += ( float )Y_[k] * digammaf( alpha );
	// the forth term
	if( k == 0 ){
		;
	} else{
		for( int y = 0; y < Y_[k+1]; y++ ){
			int y2 = y % Y_[k];
			// the first term of the inner part
			gradient_logPostAlphas += v[k-1][y2][j] * ( digammaf( ( float )n_z_[k][y][j] + alpha * v[k-1][y2][j] ) - digammaf( alpha * v[k-1][y2][j] ) );
		}
		// the second term of the inner part
		for( int y = 0; y < Y_[k]; y++ ){
			if( j == 0){
				gradient_logPostAlphas -= digammaf( alpha );
			} else {
				gradient_logPostAlphas -= digammaf( ( float )n_z_[k][y][j-1] + alpha );
			}
		}
	}

	return gradient_logPostAlphas;
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
				ofile_n << std::scientific << n_z_[k][y][j] << ' ';
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

