#include "CGS.h"

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
		r_[n] = ( float* )calloc( LW1, sizeof( float ) );
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
	pos_ = ( float* )calloc( LW1, sizeof( float ) );

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
	int W = motif_->getW();
	int K = Global::modelOrder;
	float v_diff;										// model parameter difference before and after iteration
	float** v_before;									// model parameter with highest order before iteration
	int y, j;

	// allocate memory for parameters v[y][j] with the highest order
	v_before = ( float** )calloc( Y_[K+1], sizeof( float* ) );
	for( y = 0; y < Y_[K+1]; y++ ){
		v_before[y] = ( float* )calloc( W, sizeof( float ) );
	}

	// iterate over
	while( iterate && ( CGSIterations_ < Global::maxCGSIterations ) ){

		CGSIterations_++;

		if( Global::verbose ){
			std::cout << CGSIterations_ << " iteration:\t";
		}

		// update model parameter with the highest order before each iteration
		for( y = 0; y < Y_[K+1]; y++ ){
			for( j = 0; j < W; j++ ){
				v_before[y][j] = motif_->getV()[K][y][j];
			}
		}

		// sampling z and q
		sampling_z_q();

		// update kmers counts n for the complete sequence set
		n_z_ = Global::posSequenceSet->countKmers( K, W, z_ );
		// update v_prior
		motif_->updateV( n_z_, alpha_ );

		// * optional: optimize hyper-parameter alpha
		if( !Global::noAlphaSampling )	alphaSampling();

		// check model parameter difference for convergence
		v_diff = 0.0f;
		for( y = 0; y < Y_[K+1]; y++ ){
			for( j = 0; j < W; j++ ){
				v_diff += fabsf( motif_->getV()[K][y][j] - v_before[y][j] );
			}
		}
		if( Global::verbose ){
			std::cout << "para_diff = " << v_diff << std::endl;
		}

		if( v_diff < Global::epsilon )		iterate = false;
	}

	// calculate probabilities
	motif_->calculateP();

	fprintf( stdout, "\n--- Runtime for Collapsed Gibbs sampling: %.4f seconds ---\n", ( ( float )( clock() - t0 ) ) / CLOCKS_PER_SEC );
}

void CGS::sampling_z_q(){

	int N = Global::posSequenceSet->getN();
	std::vector<Sequence*> posSeqs = Global::posSequenceSet->getSequences();
	int W = motif_->getW();
	int K = Global::modelOrder;
	int K_bg = Global::bgModelOrder;
	int n, y, i, j;
	int N_0 = 0;							// counts of sequences which do not contain motifs.
	// loop over all sequences and drop one sequence each time and update r
	for( n = 0; n < N; n++ ){
		int LW1 = posSeqs[n]->getL() - W + 1;

		// count k-mers occurring at the motif at position i except the n'th sequence
		n_z_ = Global::posSequenceSet->countKmers( K, W, z_, n );

		// updated model parameters v excluding the n'th sequence
		motif_->updateV( n_z_, alpha_ );

		// sampling equation: calculate responsibilities over all LW1 positions on n'th sequence
		std::vector<float> posterior_array;
		for( i = 0; i < LW1; i++ ){
			float posterior = 1.0f;
			for( j = 0; j < W; j++ ){
				// extract k-mers on the motif at position i over W of the n'th sequence
				y = posSeqs[n]->extractKmer( i+j, std::min( i+j, K ) );
				posterior *= ( motif_->getV()[K][y][j] / bg_->getV()[std::min( K, K_bg )][y] );
			}
			posterior_array.push_back( posterior );
		}

		// draw a new position z from discrete posterior distribution
		std::discrete_distribution<> posterior_dist( posterior_array.begin(), posterior_array.end() );
		std::default_random_engine rand;		// pick a random number
		z_[n] = posterior_dist( rand );			// draw the sample

		// checking z values from the first 10 sequences
		if( Global::verbose ){
			if( n < 10 ) std::cout << z_[n] << '\t';
		}

		if( z_[n] == 0 ) N_0++;
	}
	// sampling q_:
	// draw two random numbers Q and P from Gamma distribution
	std::gamma_distribution<> P_Gamma_dist( N_0 + 1, 1 );
	std::gamma_distribution<> Q_Gamma_dist( N - N_0 + 1, 1 );
	std::default_random_engine rand1;		// pick a random number
	std::default_random_engine rand2;		// pick another random number
	double P = P_Gamma_dist( rand1 );		// draw a sample for P
	double Q = Q_Gamma_dist( rand2 );		// draw a sample for Q

	q_ = Q / ( Q + P );						// calculate q_

	if( Global::verbose ){
		// checking z values from the first 10 sequences
		if( n < 10 ) std::cout << z_[n] << '\t';
		std::cout << N_0 << " sequences do not have motif. q = " << q_  << "\t";
	}
}

void CGS::alphaSampling(){
	// update alphas due to the learning rate eta and gradient of the log posterior of alphas

//	int K = Global::modelOrder;
//	int W = motif_->getW();

	// calcGrad_logPostAlphas();

}

float CGS::calcGrad_logPostAlphas( float alpha, int k ){
	// calculate gradient of the log posterior of alphas
	float gradient_logPostAlphas;
	int W = motif_->getW();
	float*** v = motif_->getV();

	// the first term
	gradient_logPostAlphas = ( -2.0f ) / alpha;
	// the second term
	gradient_logPostAlphas += Global::modelBeta * powf( Global::modelGamma, ( float )k ) / powf( alpha, 2.0f );
	// the third term
	gradient_logPostAlphas += ( float )Y_[k] * ( float )W * digammaf( alpha );
	// the forth term
	for( int j = 0; j < W; j++){
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

/*
	// output responsibilities r[n][i]
	std::string opath_r = opath + ".CGSposterior";
	std::ofstream ofile_r( opath_r.c_str() );
	for( int n = 0; n < Global::posSequenceSet->getN(); n++ ){
		for( int i = 0; i < posSeqs[n]->getL()-W+1; i++ ){
			ofile_r << std::scientific << r_[n][i] << ' ';
		}
		ofile_r << std::endl;
	}
*/

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

