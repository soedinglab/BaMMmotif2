/*
 * ModelLearning.cpp
 *
 *  Created on: Dec 2, 2016
 *      Author: wanwan
 */

#include "ModelLearning.h"
#include "SeqGenerator.h"

#include <cmath>								/* lgamma function */

#include <boost/math/special_functions.hpp>		/* gamma function and digamma function */
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/random.hpp>
#include <cassert>
#include <stdlib.h>

ModelLearning::ModelLearning( Motif* motif, BackgroundModel* bg, std::vector<int> folds ){

	motif_ = motif;
	bg_ = bg;

	q_ = Global::q;

	// define parameters
	int K = Global::modelOrder;
	int W = motif_->getW();
	int y, k, j, LW2;
	for( k = 0; k < Global::Yk; k++ ){
		Y_.push_back( ipow( Alphabet::getSize(), k ) );
	}

	// deep copy positive sequences
	if( folds.empty() ){										// for EM or Gibbs sampling
		folds_.resize( Global::cvFold );
		std::iota( std::begin( folds_ ), std::end( folds_ ), 0 );
		posSeqs_ = Global::posSequenceSet->getSequences();
	} else {													// for cross-validation
		folds_ = folds;
		// copy selected positive sequences due to folds
		posSeqs_.clear();
		for( size_t f_idx = 0; f_idx < folds_.size(); f_idx++ ){
			for( size_t s_idx = 0; s_idx < Global::posFoldIndices[folds_[f_idx]].size(); s_idx++ ){
				posSeqs_.push_back( Global::posSequenceSet->getSequences()[Global::posFoldIndices[folds_[f_idx]][s_idx]] );
			}
		}
	}

	// allocate memory for r_[n][i], pos_[n][i], z_[n], normFactor_[n] and initialize them
	r_ = ( float** )calloc( posSeqs_.size(), sizeof( float* ) );
	pos_ = ( float** )calloc( posSeqs_.size(), sizeof( float* ) );
	z_ = ( int* )calloc( posSeqs_.size(), sizeof( int ) );

	for( size_t n = 0; n < posSeqs_.size(); n++ ){
		LW2 = posSeqs_[n]->getL() - W + 2;
		r_[n] = ( float* )calloc( LW2, sizeof( float ) );
		pos_[n] = ( float* )calloc( LW2, sizeof( float ) );
	}

	// allocate memory for n_[k][y][j] and probs_[k][y][j] and initialize them
	n_ = ( float*** )calloc( K+1, sizeof( float** ) );
	n_z_ = ( int*** )calloc( K+1, sizeof( int** ) );
	for( k = 0; k < K+1; k++ ){
		// allocate -1 position for y
		n_[k] = ( float** )calloc( Y_[k+1]+1, sizeof( float* ) )+1;
		n_z_[k] = ( int** )calloc( Y_[k+1]+1, sizeof( int* ) )+1;
		// allocate -K positions for j
		for( y = -1; y < Y_[k+1]; y++ ){
			n_[k][y] = ( float* )calloc( W + K, sizeof( float ) )+K;
			n_z_[k][y] = ( int* )calloc( W + K, sizeof( int ) )+K;
		}
	}

	// allocate memory for alpha_[k][j] and initialize it
	alpha_ = ( double** )malloc( ( K+1 ) * sizeof( double* ) );
	m1_t_ = ( double** )calloc( K+1, sizeof( double* ) );
	m2_t_ = ( double** )calloc( K+1, sizeof( double* ) );
	for( k = 0; k < K+1; k++ ){
		alpha_[k] = ( double* )malloc( W * sizeof( double ) );
		m1_t_[k] = ( double* )calloc( W, sizeof( double ) );
		m2_t_[k] = ( double* )calloc( W, sizeof( double ) );
		for( j = 0; j < W; j++ ){
			alpha_[k][j] = Global::modelAlpha[k];
		}
	}

}

ModelLearning::~ModelLearning(){

	int K = Global::modelOrder;

	for( size_t n = 0; n < posSeqs_.size(); n++ ){
		free( r_[n] );
		free( pos_[n] );
	}
	free( r_ );
	free( pos_ );
	free( z_ );

	for( int k = 0; k < K+1; k++ ){
		for( int y = -1; y < Y_[k+1]; y++ ){
			free( n_[k][y] - K );
			free( n_z_[k][y] - K );
		}
		free( n_[k]-1 );
		free( n_z_[k]-1 );
		free( alpha_[k] );
		free( m1_t_[k] );
		free( m2_t_[k] );
	}
	free( n_ );
	free( n_z_ );
	free( alpha_ );
	free( m1_t_ );
	free( m2_t_ );

}

int ModelLearning::EM(){

	fprintf( stderr," ______\n"
					"|      |\n"
					"|  EM  |\n"
					"|______|\n\n" );

	clock_t t0 = clock();
	bool 	iterate = true;									// flag for iterating before convergence
	int 	W = motif_->getW();
	int 	K = Global::modelOrder;

	int 	y, j;
	float 	v_diff, llikelihood_prev, llikelihood_diff = 0.0f;
	float**	v_before;										// hold the parameters of the highest-order before EM

	// allocate memory for parameters v[y][j] with the highest order
	v_before = ( float** )calloc( Y_[K+1], sizeof( float* ) );
	for( y = 0; y < Y_[K+1]; y++ ){
		v_before[y] = ( float* )calloc( W, sizeof( float ) );
	}

	int EMIterations = 0;
	// iterate over
	while( iterate && ( EMIterations < Global::maxEMIterations ) ){

		EMIterations++;

		// get parameter variables with highest order before EM
		llikelihood_prev = llikelihood_;
		for( y = 0; y < Y_[K+1]; y++ ){
			for( j = 0; j < W; j++ ){
				v_before[y][j] = motif_->getV()[K][y][j];
			}
		}

		// E-step: calculate posterior
		EStep();

		// M-step: update model parameters
		MStep();

		// * optional: optimize parameter q
//		if( !Global::noQOptimization )		optimize_q();

		// check parameter difference for convergence
		v_diff = 0.0f;
		for( y = 0; y < Y_[K+1]; y++ ){
			for( j = 0; j < W; j++ ){
				v_diff += fabsf( motif_->getV()[K][y][j] - v_before[y][j] );
			}
		}

		// check the change of likelihood for convergence
		llikelihood_diff = llikelihood_ - llikelihood_prev;

		if( Global::verbose ){
			std::cout << EMIterations << " iteration:	";
			std::cout << "para_diff = " << v_diff << ",	";
			std::cout << "log likelihood = " << llikelihood_ << " 	";
			if( llikelihood_diff < 0 && EMIterations > 1 ) std::cout << " decreasing... ";
			std::cout << std::endl;
		}

		if( v_diff < Global::epsilon )					iterate = false;
		if( llikelihood_diff < 0 && EMIterations > 1 )	iterate = false;
	}

	// calculate probabilities
	motif_->calculateP();

	// free memory
	for( y = 0; y < Y_[K+1]; y++ ){
		free( v_before[y] );
	}
	free( v_before );

	fprintf( stdout, "\n--- Runtime for EM: %.4f seconds ---\n", ( ( float )( clock() - t0 ) ) / CLOCKS_PER_SEC );

    return 0;
}

void ModelLearning::EStep(){

	int K = Global::modelOrder;
	int W = motif_->getW();
	llikelihood_ = 0.0f;

	motif_->calculateLinearS( bg_->getV() );
	float** s = motif_->getS();
	int N0 = 0;
	// parallel the code
//	#pragma omp parallel for

	// calculate responsibilities r_[n][i] at position i in sequence n
	// n runs over all sequences
	for( size_t n = 0; n < posSeqs_.size(); n++ ){

		int L = posSeqs_[n]->getL();
		int LW1 = L - W + 1;
		int LW2 = L - W + 2;
		int* kmer = posSeqs_[n]->getKmer();
		float normFactor = 0.0f;

		// reset r_[n][i] and pos_[n][i]
		for( int i = 0; i < LW2; i++ ){
			r_[n][i] = 1.0f;
			pos_[n][i] = q_ / static_cast<float>( LW1 );
		}
		pos_[n][0] = 1 - q_;

		// when p(z_n > 0)
		// ij = i+j runs over all positions in sequence
		for( int ij = 0; ij < L; ij++ ){

			// extract (K+1)-mer y from positions (i-k,...,i)
			int y = ( kmer[ij] >= 0 ) ? kmer[ij] % Y_[K+1] : -1;
			// j runs over all motif positions
			for( int j = ( 0 > ( ij-L+W ) ? 0 : ij-L+W ); j < ( W < (ij+1) ? W : ij+1 ); j++ ){
				r_[n][LW1-ij+j] *= s[y][j];
			}
		}

		// calculate complete responsibilities and sum them up
		for( int i = 0; i < LW2; i++ ){
			r_[n][i] *= pos_[n][i];
			normFactor += r_[n][i];
		}

		// normalize responsibilities
		for( int i = 0; i < LW2; i++ ){
			r_[n][i] /= normFactor;
		}

		if( r_[n][0] > 0.7 ) N0++;
		// calculate log likelihood over all sequences
		llikelihood_ += logf( normFactor );
	}

//	std::cout << "N0=" << N0 << std::endl;
}

void ModelLearning::MStep(){

	int K = Global::modelOrder;
	int W = motif_->getW();

	// reset the fractional counts n
	for( int k = 0; k < K+1; k++ ){
		for( int y = 0; y < Y_[k+1]; y++ ){
			for( int j = 0; j < W; j++ ){
				n_[k][y][j] = 0.0f;
			}
		}
	}

	// compute fractional occurrence counts for the highest order K
	// n runs over all sequences
	for( size_t n = 0; n < posSeqs_.size(); n++ ){
		int L = posSeqs_[n]->getL();
		int* kmer = posSeqs_[n]->getKmer();

		// ij = i+j runs over all positions i on sequence n
		for( int ij = 0; ij < L; ij++ ){

			int y = ( kmer[ij] >= 0 ) ? kmer[ij] % Y_[K+1] : -1;

			for( int j = ( 0 > ( ij-L+W ) ? 0 : ij-L+W ); j < ( W < (ij+1) ? W : ij+1 ); j++ ){

				n_[K][y][j] += r_[n][L-W+1-ij+j];

			}
		}
	}

	// compute fractional occurrence counts from higher to lower order
	// k runs over all orders
	for( int k = K; k > 0; k-- ){

		for( int y = 0; y < Y_[k+1]; y++ ){

			int y2 = y % Y_[k];

			for( int j = 0; j < W; j++ ){

				n_[k-1][y2][j] += n_[k][y][j];

			}
		}
	}

	// update model parameters v[k][y][j]
	motif_->updateV( n_, alpha_, K );
}

void ModelLearning::optimize_q(){

	// optimize hyper-parameter q
	// motif.updateV()
}

void ModelLearning::GibbsSampling(){

	fprintf(stderr, " ___________________________\n"
					"|                           |\n"
					"|  Collapsed Gibbs sampler  |\n"
					"|___________________________|\n\n" );

	clock_t t0 = clock();
	int iteration = 0;

	int K = Global::modelOrder;
	int W = motif_->getW();
	float eta = Global::eta;						// learning rate for alpha learning
	int k, y, j, n, i;

	// initialize z for all the sequences
	if( !Global::noInitialZ ){

		// E-step: calculate posterior
		EStep();

		// extract initial z from the indices of the largest responsibilities
		for( n = 0; n < ( int )posSeqs_.size(); n++ ){
			int LW1 = posSeqs_[n]->getL() - W + 1;
			float maxR = r_[n][0];
			int maxIdx = 0;
			for( i = 1; i <= LW1; i++ ){
				if( r_[n][i] > maxR ){
					maxR = r_[n][i];
					maxIdx = LW1+1-i;
				}
			}
			z_[n] = maxIdx;
		}

	} else {
		// initialize z with a random number
		for( n = 0; n < ( int )posSeqs_.size(); n++ ){
			int LW2 = posSeqs_[n]->getL() - motif_->getW() + 2;
			z_[n] = rand() % LW2;
		}
	}

	// count the k-mers
	// 1. reset n_z_[k][y][j] to 0
	for( k = 0; k < K+1; k++ ){
		for( y = 0; y < Y_[k+1]; y++ ){
			for( j = -K; j < W; j++ ){
				n_z_[k][y][j] = 0;
			}
		}
	}

	// 2. count k-mers for the highest order K
	for( n = 0; n < ( int )posSeqs_.size(); n++ ){
		if( z_[n] > 0 ){
			int* kmer = posSeqs_[n]->getKmer();
			for( j = ( z_[n] <= K ) ? 1-z_[n] : -K; j < W; j++ ){
				y = ( kmer[z_[n]-1+j] >= 0 ) ? kmer[z_[n]-1+j] % Y_[K+1] : -1;
				n_z_[K][y][j]++;
			}
		}
	}

	// compute k-mer counts for all the lower orders
	for( k = K; k > 0; k-- ){
		for( y = 0; y < Y_[k+1]; y++ ){
			int y2 = y % Y_[k];
			for( j = -K; j < W; j++ ){
				n_z_[k-1][y2][j] += n_z_[k][y][j];
			}
		}
	}

	// Gibbs sampling position z and fraction q
	for( int iter = 0; iter < 10; iter++ ){
		Gibbs_sample_z_q();
	}

	// vector to store the last a few alphas for sampling methods
	std::vector<std::vector<double>> alpha_avg( K+1 );
	for( k = 0; k < K+1; k++ ){
		alpha_avg[k].resize( W );
		for( j = 0; j < W; j++ ){
			alpha_avg[k][j] = 0.0;
		}
	}

	// todo: only for writing out the log posterior of alphas
	std::string opath = std::string( Global::outputDirectory ) + "/Kj.lposA";
	std::ofstream ofile( opath );


	double llikelihood = calc_llikelihood_alphas( alpha_, K );
	// iterate over
	while( iteration < Global::maxCGSIterations ){

		iteration++;

//		float llikelihood_prev;
//		llikelihood_prev = llikelihood_;

		// Gibbs sampling position z and fraction q
		if( !Global::noZSampling )	Gibbs_sample_z_q();

//		std::cout << "diff_llikelihood=" << llikelihood_ - llikelihood_prev << std::endl;

		// update alphas by stochastic optimization
		if( !Global::noAlphaOptimization ){

			stochastic_optimize_alphas( K, W, eta, iteration );

/*			for( j = 0; j < W; j++ ){

				ofile << std::scientific << std::setprecision( 8 ) << calc_gradient_alphas( alpha_, K, j ) << '\t' << alpha_[K][j] << '\t';

			}
			ofile << std::endl;*/

		} else if( Global::GibbsMHalphas ){

			GibbsMH_sample_alphas( iteration );

			if( iteration > Global::maxCGSIterations - 10 ){
				for( k = 0; k < K+1; k++ ){
					for( j = 0; j < W; j++ ){
						alpha_avg[k][j] += alpha_[k][j];
					}
				}
			}

/*			if( iteration > 35 ){

				for( j = 0; j < W; j++ ){

					ofile << std::scientific << std::setprecision( 8 ) << calc_logCondProb_a( iteration, log( alpha_[K][j] ), K, j )
							<< '\t' << log( alpha_[K][j] )<< '\t';
				}

				ofile << std::endl;
			}*/

		} else if( Global::dissampleAlphas ){

			discrete_sample_alphas( iteration );

			if( iteration > Global::maxCGSIterations - 10 ){
				for( k = 0; k < K+1; k++ ){
					for( j = 0; j < W; j++ ){
						alpha_avg[k][j] += alpha_[k][j];
					}
				}
			}
/*			// todo:
			for( j = 0; j < W; j++ ){

				ofile << std::scientific << std::setprecision( 8 ) << calc_logCondProb_a( iteration, log( alpha_[K][j] ), K, j )
						<< '\t' << log( alpha_[K][j] )<< '\t';
			}
			ofile << std::endl;*/

		} else if( Global::debugAlphas ){

			// todo: only for writing out the log posterior of alphas
			for( j = 0; j < W; j++ ){
				ofile << std::scientific << std::setprecision( 8 ) << calc_logCondProb_a( iteration, ( double )iteration / 20.0, K, j )
						<< '\t' << ( double )iteration / 20.0 << '\t';
			}
			ofile << std::endl;
		}

		// todo: check the log likelihood
		std::cout << calc_llikelihood_alphas( alpha_, K ) - llikelihood << std::endl;
		llikelihood = calc_llikelihood_alphas( alpha_, K );
	}

	// obtaining a motif model:
	if( Global::GibbsMHalphas || Global::dissampleAlphas ){
		// average alphas over the last few steps for GibbsMH
		for( k = 0; k < K+1; k++ ){
			for( j = 0; j < W; j++ ){
				alpha_[k][j] = alpha_avg[k][j] / 10.0;
			}
		}
	}
	// update model parameter v
	motif_->updateVz_n( n_z_, alpha_, K );

	// run five steps of EM to optimize the final model with
	// the optimum model parameters v's and the fixed alphas
	for( size_t step = 0; step < 5; step++ ){

		// E-step: calculate posterior
		EStep();

		// M-step: update model parameters
		MStep();
	}

	// calculate probabilities
	motif_->calculateP();

	// update the global parameter q
	Global::q = q_;

	fprintf( stdout, "\n--- Runtime for Collapsed Gibbs sampling: %.4f seconds ---\n", ( ( float )( clock() - t0 ) ) / CLOCKS_PER_SEC );
}

void ModelLearning::Gibbs_sample_z_q(){

	int N = static_cast<int>( posSeqs_.size() );
	int W = motif_->getW();
	int K = Global::modelOrder;
	int K_bg = ( Global::bgModelOrder < K ) ? Global::bgModelOrder : K;
	int k, y, j, i, n;
	int N_0 = 0;
	float N_1 = static_cast<float>( posSeqs_.size() ) * q_;

	llikelihood_ = 0.0f;

	float*** v = motif_->getV();
	float** v_bg = bg_->getV();

	// updated model parameters v excluding the n'th sequence
	motif_->updateVz_n( n_z_, alpha_, K );
	// compute log odd scores s[y][j], log likelihoods of the highest order K
	motif_->calculateLinearS( v_bg );

	float** s = motif_->getS();

	// sampling z:
	// loop over all sequences and drop one sequence each time and update r
	for( n = 0; n < N; n++ ){

		int L = posSeqs_[n]->getL();
		int LW1 = L - W + 1;
		int* kmer = posSeqs_[n]->getKmer();

		// count k-mers at position z_i+j except the n'th sequence
		// remove the k-mer counts from the sequence with the current z
		// and re-calculate the log odds scores in linear space
		/*
		 * -------------- faster version of removing k-mer ------------------
		 */

		if( z_[n] > 0 ){

			for( k = 0; k < K+1; k++ ){

				for( j = ( z_[n] <= K ) ? 1-z_[n] : -K; j < W; j++ ){

					y = ( kmer[z_[n]-1+j] >= 0 ) ? kmer[z_[n]-1+j] % Y_[k+1] : -1;

					n_z_[k][y][j]--;

					if( j >= 0 && y >= 0 ){
						if( k == 0 ){
							v[k][y][j]= ( static_cast<float>( n_z_[k][y][j] ) + static_cast<float>( alpha_[k][j] ) * v_bg[k][y] )
										/ ( N_1 + static_cast<float>( alpha_[k][j] ) );
						} else {

							int y2 = y % Y_[k];
							int yk = y / Y_[1];

							v[k][y][j] = ( static_cast<float>( n_z_[k][y][j] ) + static_cast<float>( alpha_[k][j] ) * v[k-1][y2][j] )
										/ ( static_cast<float>( n_z_[k-1][yk][j-1] ) + static_cast<float>( alpha_[k][j] ) );
						}

						if( k == K ){

							int y_bg = y % Y_[K_bg+1];

							s[y][j] = v[K][y][j] / v_bg[K_bg][y_bg];

						}
					}
				}
			}
		}

		/*
		 * -------------- slower version of removing k-mer -----------------------
		 */
/*

		if( z_[n] > 0 ){
			for( j = ( z_[n] <= K ) ? 1-z_[n] : -K; j < W; j++ ){
				for( k = 0; k < K+1; k++ ){
					y = ( kmer[z_[n]-1+j] >= 0 ) ? kmer[z_[n]-1+j] % Y_[k+1] : -1;
					n_z_[k][y][j]--;
				}
			}
		}

		// updated model parameters v excluding the n'th sequence
		motif_->updateVz_n( n_z_, alpha_, K );
		// compute log odd scores s[y][j], log likelihoods of the highest order K
		motif_->calculateLinearS( bg_->getV() );

		s = motif_->getS();
*/

		/*
		 * ------- sampling equation -------
		 */

		// calculate responsibilities and positional priors
		// over all LW1 positions on n'th sequence:
		float normFactor = 0.0f;
		pos_[n][0] = 1.0f - q_;
		r_[n][0] = pos_[n][0];
		normFactor += r_[n][0];
		float pos_i = q_ / static_cast<float>( LW1 );
		for( i = 1; i <= LW1; i++ ){
			pos_[n][i] = pos_i;
			r_[n][i] = 1.0f;
		}

		// todo: could be parallelized by extracting 8 sequences at once
		// ij = i+j runs over all positions in sequence
		for( int ij = 0; ij < L; ij++ ){

			// extract (K+1)-mer y from positions (i-k,...,i)
			y = ( kmer[ij] >= 0 ) ? kmer[ij] % Y_[K+1] : -1;

			// j runs over all motif positions
			for( j = ( ( ij-L+W ) < 0 ? 0 : ij-L+W ); j < ( W < (ij+1) ? W : ij+1 ); j++ ){
				r_[n][LW1-ij+j] *= s[y][j];
			}
		}
		for( i = 1; i <= LW1; i++ ){
			r_[n][i] *= pos_[n][i];
			normFactor += r_[n][i];
		}

		// calculate log likelihood of sequences given the motif positions z and model parameter v
		llikelihood_ += logf( normFactor );

		// normalize responsibilities and append them to an array
		std::vector<float> posteriors;
		r_[n][0] /= normFactor;
		posteriors.push_back( r_[n][0] );

		for( i = LW1; i >= 1; i-- ){
			r_[n][i] /= normFactor;
			posteriors.push_back( r_[n][i] );
		}

		// draw a new position z from the discrete distribution of posterior
		std::discrete_distribution<> posterior_dist( posteriors.begin(), posteriors.end() );
		z_[n] = posterior_dist( Global::rngx );

		if( z_[n] == 0 ){
			// count sequences which do not contain motifs.
			N_0++;

		} else {
			// add the k-mer counts from the current sequence with the updated z
			for( j = ( z_[n] <= K ) ? 1-z_[n] : -K; j < W; j++ ){
				for( k = 0; k < K+1; k++ ){
					y = ( kmer[z_[n]-1+j] >= 0 ) ? kmer[z_[n]-1+j] % Y_[k+1] : -1;
					n_z_[k][y][j]++;
				}
			}
		}

	}

	// sampling q:
	if( !Global::noQSampling ){
		boost::math::beta_distribution<> q_beta_dist( N - N_0 + 1, N_0 + 1);
		q_ = ( float )quantile( q_beta_dist, ( float )rand() / ( float )RAND_MAX );
	}

//	std::cout << "N=" << N - N_0 << ",\t q=" << ( float )( N - N_0 ) / ( float )N  << ",\t";

}

void ModelLearning::stochastic_optimize_alphas( int K, int W, float eta, int iter ){
	// update alphas using stochastic optimization algorithm ADAM (DP Kingma & JL Ba 2015)

	int k, j;
	double beta1 = 0.9;		// exponential decay rate for the moment estimates
	double beta2 = 0.99;	// exponential decay rate for the moment estimates
	double epsilon = 1e-8;	// cutoff
	double gradient;		// gradient of log posterior of alpha
	double m1;				// first moment vector (the mean)
	double m2;				// second moment vector (the uncentered variance)

	double t = static_cast<double>( iter );

	for( k = 0; k < K+1; k++ ){

		for( j = 0; j < W; j++ ){

			// re-parameterise alpha on log scale: alpha = e^a
			double a  = log( alpha_[k][j] );

			// get gradients w.r.t. stochastic objective at timestep t
			gradient = alpha_[k][j] * calc_gradient_alphas( alpha_, k, j );

			// update biased first moment estimate
			m1_t_[k][j] = beta1 * m1_t_[k][j] + ( 1 - beta1 ) * gradient;

			// update biased second raw moment estimate
			m2_t_[k][j] = beta2 * m2_t_[k][j] + ( 1 - beta2 ) * gradient * gradient;

			// compute bias-corrected first moment estimate
			m1 = m1_t_[k][j] / ( 1 - pow( beta1, t ) );

			// compute bias-corrected second raw moment estimate
			m2 = m2_t_[k][j] / ( 1 - pow( beta2, t ) );

			// update parameter a due to alphas
			// Note: here change the sign in front of eta from '-' to '+'
			a += eta * m1 / ( ( sqrt( m2 ) + epsilon ) * sqrt( t ) );

			alpha_[k][j] = exp( a );
		}

	}

}

void ModelLearning::GibbsMH_sample_alphas( int iter ){
	// Gibbs sampling alphas in exponential space with Metropolis-Hastings algorithm

	int K = Global::modelOrder;
	int W = motif_->getW();

	for( int k = 0; k < K+1; k++ ){

		for( int j = 0; j < W; j++ ){

			// draw 10 times in a row and take record of the last accepted sample
//			for( int step = 0; step < 10; step++ ){

				// Metropolis-Hasting sheme
				double a_prev = log( alpha_[k][j] );

				double lprob_a_prev = calc_logCondProb_a( iter, a_prev, k, j );

				// draw a new 'a' from the distribution of N(a, 1)
//				std::normal_distribution<float> norm_dist( a_prev, 1.0f );
//				std::normal_distribution<float> norm_dist( a_prev, 1.0f / sqrtf( ( float )( k+1 ) ) );
				std::normal_distribution<double> norm_dist( a_prev, 1.0 / ( double )( k+1 ) );

				double a_new = norm_dist( Global::rngx );

				double lprob_a_new = calc_logCondProb_a( iter, a_new, k, j );
				double accept_ratio;
				double uni_random;
				if( lprob_a_new < lprob_a_prev ){
					// calculate the acceptance ratio
					accept_ratio = exp( lprob_a_new - lprob_a_prev );

					// draw a random number uniformly between 0 and 1
					std::uniform_real_distribution<double> uniform_dist( 0.0, 1.0 );
					uni_random = uniform_dist( Global::rngx );

					// accept the trial sample if the ratio is not smaller than a random number between (0,1)
					if( accept_ratio >= uni_random ){

						alpha_[k][j] = exp( a_new );

					}

				} else {
					// accept the trial sample
					alpha_[k][j] = exp( a_new );

				}

/*				if( k == K && j == 4 ) std::cout << iter << ": j=5, "
										<< std::setprecision( 8 )
										<< "ao=" << a_prev
										<< ",\t po=" << lprob_a_prev
										<< ",\t an=" << a_new
										<< ",\t pn=" << lprob_a_new
										<< ",\t Racc=" << accept_ratio
										<< ",\t Rrej=" << uni_random
										<< std::endl;*/
			}
//		}
	}
}

void ModelLearning::discrete_sample_alphas( int iter ){
	// sample an alpha from the discrete distribution of its log posterior

	int K = Global::modelOrder;
	int W = motif_->getW();

	std::string opath = std::string( Global::outputDirectory ) + '/'
						+ Global::posSequenceBasename + ".lpobKj";
	std::ofstream ofile( opath.c_str() );

	for( int k = 0; k < K+1; k++ ){

		for( int j = 0; j < W; j++ ){

			std::vector<double> condProb;

			double condProb_new;

			double base = calc_logCondProb_a( iter, 0.0, k, j );

			for( int it = 0; it < 100; it++ ){

				condProb_new = exp( calc_logCondProb_a( iter, static_cast<double>( it ) / 10.0, k, j ) - base );

				condProb.push_back( condProb_new );

			}

			std::discrete_distribution<> posterior_dist( condProb.begin(), condProb.end() );

			alpha_[k][j] = exp( posterior_dist( Global::rngx ) / 10.0 );
		}
	}
}

double ModelLearning::calc_gradient_alphas( double** alpha, int k, int j ){
	// calculate partial gradient of the log posterior of alphas due to equation 47 in the theory
	// Note that j >= k

	double gradient = 0.0;
	float*** v = motif_->getV();
	float** v_bg = bg_->getV();
	int N = static_cast<int>( posSeqs_.size() ) - 1;
	int y, y2;


	// the first term of equation 47
	gradient -= 2.0 / alpha[k][j];

	// the second term of equation 47
	gradient += Global::modelBeta * pow( Global::modelGamma, static_cast<double>( k ) ) / pow( alpha[k][j], 2.0 );

	// the third term of equation 47
	gradient += static_cast<double>( ipow( Y_[1], k ) ) * boost::math::digamma( alpha[k][j] );

	// the forth term of equation 47
	for( y = 0; y < Y_[k+1]; y++ ){

		if( k == 0 ){
			// the first part
			gradient += v_bg[k][y] * boost::math::digamma( static_cast<double>( n_z_[k][y][j] ) + alpha[k][j] * v_bg[k][y] );

			// the second part
			gradient -= v_bg[k][y] * boost::math::digamma( alpha[k][j] * v_bg[k][y] );

		} else {

			y2 = y % Y_[k];

			// the first part
			gradient += v[k-1][y2][j] * boost::math::digamma( static_cast<double>( n_z_[k][y][j] ) + alpha[k][j] * v[k-1][y2][j] );

			// the second part
			gradient -= v[k-1][y2][j] * boost::math::digamma( alpha[k][j] * v[k-1][y2][j] );

		}

	}

	// the last term of equation 47
	for( y = 0; y < Y_[k]; y++ ){

		if( k == 0 ){

			gradient -= boost::math::digamma( static_cast<double>( N ) + alpha[k][j] );

		} else if( j == 0 ){
			int sum = 0;
			for( int a = 0; a < Y_[1]; a++ ){
				int ya = y * Y_[1] + a;
				sum += n_z_[k][ya][j];
			}
			gradient -= boost::math::digamma( static_cast<double>( sum ) + alpha[k][j] );

		} else {

			gradient -= boost::math::digamma( static_cast<double>( n_z_[k-1][y][j-1] ) + alpha[k][j] );

		}
	}
	return gradient;
}

double ModelLearning::calc_logCondProb_a( int iteration, double a, int k, int j ){
	// calculate partial log conditional probabilities of a's due to equation 50 in the theory

	double logCondProbA = 0.0f;
	float*** v = motif_->getV();
	float** v_bg = bg_->getV();
	int N = ( int )posSeqs_.size() - 1;
	int y, y2, ya;

	bool be_printed = false;

/*	if( iteration == 3 && j == 12 && k == Global::modelOrder ){
		be_printed = true;
	}*/

	// get alpha by alpha = e^a
	double alpha = exp( a );

	// the first term of equation 50
	logCondProbA -= a;

	if( be_printed ) std::cout << std::endl << "at iteration " << iteration+1 << " when a=" << a << std::endl
			<< -a << '\t' << logCondProbA << "\t(1)" << std::endl;

	// the second term of equation 50
	logCondProbA -= Global::modelBeta * pow( Global::modelGamma, ( double )k ) / alpha;

	if( be_printed ) std::cout << -Global::modelBeta * pow( Global::modelGamma, ( double )k ) / alpha
								<< '\t' << logCondProbA << "\t(2)" << std::endl;

	if( k == 0 ){

		for( y = 0; y < Y_[k]; y++ ){

			// the third term of equation 50
			logCondProbA += boost::math::lgamma( alpha );

			// the first part of the forth term of equation 50
			logCondProbA += boost::math::lgamma( static_cast<double>( n_z_[k][y][j] ) + alpha * v_bg[k][y] );

			// the second part of the forth term of equation 50
			logCondProbA -= boost::math::lgamma( alpha * v_bg[k][y] );

		}

		// the fifth term of equation 50
		logCondProbA -= boost::math::lgamma( static_cast<double>( N ) + alpha );

	} else {

		// the third, forth and fifth terms of equation 50
		for( y = 0; y < Y_[k]; y++ ){

			// the third term of equation 50
			logCondProbA += boost::math::lgamma( alpha );

			if( be_printed ) std::cout << boost::math::lgamma( alpha ) << '\t' << logCondProbA << "\t(3+" << y <<") : lgamma(" << alpha << ")"<< std::endl;


			// the forth term of equation 50
			for( int A = 0; A < Y_[1]; A++ ){

				ya = y * Y_[1] + A;

				y2 = ya % Y_[k];

				if( n_z_[k][ya][j] > 0 ){

					// the first part of the forth term
					logCondProbA += boost::math::lgamma( static_cast<double>( n_z_[k][ya][j] ) + alpha * v[k-1][y2][j] );

					// the second part of the forth term
					logCondProbA -= boost::math::lgamma( alpha * v[k-1][y2][j] );
				}

				if( be_printed ) std::cout << boost::math::lgamma( static_cast<double>( n_z_[k][ya][j] ) + alpha * v[k-1][y2][j] )
											- boost::math::lgamma( alpha * v[k-1][y2][j] )
											<<'\t' << logCondProbA
											<< "\t(4+" << ya << ") : lgmma(" << ( double )n_z_[k][ya][j] << "+"
											<< alpha << "*" << v[k-1][y2][j] << ")-lgamma (" << alpha << "*"
											<< v[k-1][y2][j] <<")"<< std::endl;
			}

			if( be_printed )  std::cout << std::endl;

			// the fifth term
			logCondProbA -= boost::math::lgamma( static_cast<double>( n_z_[k-1][y][j-1] ) + alpha );

			if( be_printed ) std::cout << -boost::math::lgamma( static_cast<double>( n_z_[k-1][y][j-1] ) + alpha ) << '\t' << logCondProbA
											<< "\t(5) : " << "-lgamma("<<  static_cast<double>( n_z_[k-1][y][j-1] )
											<< "+" << alpha <<")" << std::endl;

			// !!! important: correction for the occasions when zero or one k-mer is present
/*			if( n_z_[k-1][y][j-1] <= 1 ){

				logCondProbA = -a - Global::modelBeta * pow( Global::modelGamma, static_cast<double>( k ) ) / alpha;

				if( be_printed ) std::cout << "* 6 : " << logCondProbA << " since n_[k-1][y][j-1]<=1" << std::endl;

			}

			if( be_printed ) std::cout << std::endl;*/
		}

	}

	return logCondProbA;
}

double ModelLearning::calc_prior_alphas( double** alpha, int k ){
	// calculate partial log conditional probabilities of alphas due to equation 46 in the theory

	double logPrior = 0.0;
	int W = motif_->getW();

	for( int j = 0; j < W; j++ ){

		// the first term of equation 46
		logPrior -= 2.0 * log( alpha[k][j] );

		// the second term of equation 46
		logPrior -= Global::modelBeta * pow( Global::modelGamma, ( double )k ) / alpha[k][j];

	}

	return logPrior;
}

double ModelLearning::calc_llikelihood_alphas( double** alpha, int k ){
	// calculate partial log likelihood of alphas due to equation 46 in the theory

	double logLikelihood = 0.0;
	float*** v = motif_->getV();
	float** v_bg = bg_->getV();
	double N = static_cast<double>( posSeqs_.size() ) - 1.0;
	int W = motif_->getW();
	int y, y2, ya, j;



	if( k == 0 ){
		for( j = 0; j < W; j++ ){
			for( y = 0; y < Y_[k]; y++ ){

				// the third term of equation 50
				logLikelihood += boost::math::lgamma( alpha[k][j] );

				// the first part of the forth term of equation 50
				logLikelihood += boost::math::lgamma( static_cast<double>( n_z_[k][y][j] ) + alpha[k][j] * v_bg[k][y] );

				// the second part of the forth term of equation 50
				logLikelihood -= boost::math::lgamma( alpha[k][j] * v_bg[k][y] );

			}

			// the fifth term of equation 50
			logLikelihood -= boost::math::lgamma( static_cast<double>( N ) + alpha[k][j] );
		}

	} else {
		for( j = 0; j < W; j++ ){
			// the third, forth and fifth terms of equation 50
			for( y = 0; y < Y_[k]; y++ ){

				// the third term of equation 50
				logLikelihood += boost::math::lgamma( alpha[k][j] );

				// the forth term of equation 50
				for( int A = 0; A < Y_[1]; A++ ){

					ya = y * Y_[1] + A;

					y2 = ya % Y_[k];

					if( n_z_[k][ya][j] > 0 ){

						// the first part of the forth term
						logLikelihood += boost::math::lgamma( static_cast<double>( n_z_[k][ya][j] ) + alpha[k][j] * v[k-1][y2][j] );

						// the second part of the forth term
						logLikelihood -= boost::math::lgamma( alpha[k][j] * v[k-1][y2][j] );
					}

				}

				// the fifth term
				logLikelihood -= boost::math::lgamma( static_cast<double>( n_z_[k-1][y][j-1] ) + alpha[k][j] );

			}
		}
	}

	return logLikelihood;
}

Motif* ModelLearning::getMotif(){
	return motif_;
}

void ModelLearning::print(){

}

void ModelLearning::write( int N ){

	/**
	 * 	 * save EM parameters in four flat files:
	 * (1) posSequenceBasename.EMcount:			refined fractional counts of (k+1)-mers
	 * (2) posSequenceBasename.EMweight: 		responsibilities, posterior distributions
	 * (3) posSequenceBasename.EMalpha:			optimized hyper-parameter alphas
	 * (4) posSequenceBasename.EMposition:		position of motif(s) on each sequence
	 * or
	 * 	 * save CGS parameters in four flat files:
	 * (1) posSequenceBasename.CGScount:		refined integral counts of (k+1)-mers
	 * (2) posSequenceBasename.CGSweight:		responsibilities, posterior distributions
	 * (3) posSequenceBasename.CGSalpha:		optimized hyper-parameter alphas
	 * (4) posSequenceBasename.CGSposition:		position of motif on each sequence
	 */

	int k, y, j, i;
	int W = motif_->getW();
	int K = Global::modelOrder;

	std::string opath = std::string( Global::outputDirectory ) + '/'
						+ Global::posSequenceBasename + "_motif_" + std::to_string( N+1 );

	// output (k+1)-mer counts n[k][y][j]
	std::string opath_n = opath + ".kmecounts";
	std::ofstream ofile_n( opath_n.c_str() );
	if( Global::EM ){
		for( j = 0; j < W; j++ ){
			for( k = 0; k < K+1; k++ ){
				for( y = 0; y < Y_[k+1]; y++ ){
					ofile_n << ( int )n_[k][y][j] << '\t';
				}
				ofile_n << std::endl;
			}
			ofile_n << std::endl;
		}
	} else if( Global::CGS ){
		for( j = 0; j < W; j++ ){
			for( k = 0; k < K+1; k++ ){
				for( y = 0; y < Y_[k+1]; y++ ){
					ofile_n << n_z_[k][y][j] << '\t';
				}
				ofile_n << std::endl;
			}
			ofile_n << std::endl;
		}
	}

	// output position(s) of motif(s): pos_[n][i]
	std::string opath_pos = opath + ".motifpos";
	std::ofstream ofile_pos( opath_pos.c_str() );

	ofile_pos << "seq" << '\t' << "position" << std::endl;

	if( Global::EM ){
		float cutoff = 0.3f;	// threshold for having a motif on the sequence in terms of responsibilities
		for( size_t n = 0; n < posSeqs_.size(); n++ ){

			ofile_pos << n+1 << '\t';

			int LW1 = posSeqs_[n]->getL() - W + 1;

			for( i = LW1; i > 0; i-- ){
				if( r_[n][i] >= cutoff ){
					ofile_pos << LW1-i << '\t';
				}
			}

			ofile_pos << std::endl;
		}
	} else if( Global::CGS ){
		for( size_t n = 0; n < posSeqs_.size(); n++ ){
			ofile_pos << n+1 << '\t' << z_[n]+1 <<'\t' << z_[n]+W << std::endl;
		}
	}


	// output parameter alphas alpha[k][j]
	std::string opath_alpha = opath + ".alphas";
	std::string opath_a = opath + ".logalphas";
	std::ofstream ofile_alpha( opath_alpha.c_str() );
	std::ofstream ofile_a( opath_a.c_str() );
	for( k = 0; k < K+1; k++ ){
		ofile_alpha << "> k=" << k << std::endl;
		ofile_a << "> k=" << k << std::endl;
		for( j = 0; j < W; j++ ){
			ofile_alpha << std::setprecision( 3 ) << alpha_[k][j] << '\t';
			ofile_a << std::setprecision( 3 ) << log( alpha_[k][j] ) << '\t';
		}
		ofile_alpha << std::endl;
		ofile_a << std::endl;
	}

/*	// output responsibilities r[n][i]
	std::string opath_r = opath + ".weights";
	std::ofstream ofile_r( opath_r.c_str() );
	for( size_t n = 0; n < posSeqs_.size(); n++ ){
		ofile_r << std::scientific << std::setprecision( 2 ) << r_[n][0] << ' ';
		int LW1 = posSeqs_[n]->getL() - W + 1;
		for( i = LW1; i > 0; i-- ){
			ofile_r << std::setprecision( 2 ) << r_[n][i] << ' ';
		}
		ofile_r << std::endl;
	}*/


}
