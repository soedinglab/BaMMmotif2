/*
 * ModelLearning.cpp
 *
 *  Created on: Dec 2, 2016
 *      Author: wanwan
 */

#include "ModelLearning.h"
#include "SeqGenerator.h"

#include <random>
#include <cmath>

#include <boost/math/special_functions.hpp>		/* gamma function and digamma function */
#include <boost/math/special_functions/digamma.hpp>
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

	// allocate memory for s_[y][j] and initialize it
	s_ = ( float** )calloc( Y_[K+1], sizeof( float* ) );
	for( y = 0; y < Y_[K+1]; y++ ){
		s_[y] = ( float* )calloc( W, sizeof( float ) );
	}

	// allocate memory for r_[n][i], pos_[n][i], z_[n] and initialize them
	r_ = ( float** )calloc( posSeqs_.size(), sizeof( float* ) );
	pos_ = ( float** )calloc( posSeqs_.size(), sizeof( float* ) );
	z_ = ( int* )calloc( posSeqs_.size(), sizeof( int ) );
	for( size_t n = 0; n < posSeqs_.size(); n++ ){
		LW2 = posSeqs_[n]->getL() - W + 2;
		r_[n] = ( float* )calloc( LW2, sizeof( float ) );
		pos_[n] = ( float* )calloc( LW2, sizeof( float ) );
//		z_[n] = rand() % LW2;
		z_[n] = 103; //todo: fix z for testing
	}

	// allocate memory for n_[k][y][j] and probs_[k][y][j] and initialize them
	n_ = ( float*** )calloc( K+1, sizeof( float** ) );
	n_z_ = ( int*** )calloc( K+1, sizeof( int** ) );
	for( k = 0; k < K+1; k++ ){
		n_[k] = ( float** )calloc( Y_[k+1], sizeof( float* ) );
		n_z_[k] = ( int** )calloc( Y_[k+1], sizeof( int* ) );
		for( y = 0; y < Y_[k+1]; y++ ){
			n_[k][y] = ( float* )calloc( W, sizeof( float ) );
			n_z_[k][y] = ( int* )calloc( W, sizeof( int ) );
		}
	}

	// allocate memory for alpha_[k][j] and initialize it
	// allocate memory for m1_t_ and m2_t_
	alpha_ = ( float** )malloc( ( K+1 ) * sizeof( float* ) );
	m1_t_ = ( float** )calloc( K+1, sizeof( float* ) );
	m2_t_ = ( float** )calloc( K+1, sizeof( float* ) );
	for( k = 0; k < K+1; k++ ){
		alpha_[k] = ( float* )malloc( W * sizeof( float ) );
		m1_t_[k] = ( float* )calloc( W, sizeof( float ) );
		m2_t_[k] = ( float* )calloc( W, sizeof( float ) );
		for( j = 0; j < W; j++ ){
			alpha_[k][j] = Global::modelAlpha[k];
		}
	}

}

ModelLearning::~ModelLearning(){

	for( int y = 0; y < Y_[Global::modelOrder+1]; y++ ){
		free( s_[y] );
	}
	free( s_ );

	for( size_t n = 0; n < posSeqs_.size(); n++ ){
		free( r_[n] );
		free( pos_[n] );
	}
	free( r_ );
	free( pos_ );
	free( z_ );

	for( int k = 0; k < Global::modelOrder+1; k++ ){
		for( int y = 0; y < Y_[k+1]; y++ ){
			free( n_[k][y] );
			free( n_z_[k][y] );
		}
		free( n_[k] );
		free( n_z_[k] );
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
	bool iterate = true;									// flag for iterating before convergence
	int W = motif_->getW();
	int K = Global::modelOrder;
	int K_bg = ( Global::bgModelOrder < K ) ? Global::bgModelOrder : K;

	int y, y_bg, j;
	float v_diff, llikelihood_prev, llikelihood_diff = 0.0f;
	float** v_before;										// hold the parameters of the highest-order before EM

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

		// compute log odd scores s[y][j], log likelihoods of the highest order K
		for( y = 0; y < Y_[K+1]; y++ ){
			for( j = 0; j < W; j++ ){
				y_bg = y % Y_[K_bg+1];
				s_[y][j] = motif_->getV()[K][y][j] / bg_->getV()[K_bg][y_bg];
			}
		}

		// E-step: calculate posterior
		EStep();

		// M-step: update model parameters
		MStep();

		// * optional: optimize parameter q
		if( !Global::noQOptimization )		optimize_q();

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
			if( llikelihood_diff < 0 && EMIterations > 1) std::cout << " decreasing... ";
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
	llikelihood_ = 0.0f;										// reset log likelihood

	// calculate responsibilities r_[n][i] at position i in sequence n
	// n runs over all sequences
	for( size_t n = 0; n < posSeqs_.size(); n++ ){

		int L = posSeqs_[n]->getL();
		int LW1 = L - W + 1;
		int LW2 = L - W + 2;
		float normFactor = 0.0f;								// reset normalization factor

		// reset r_[n][i] and pos_[n][i]
		for( int i = 0; i < LW2; i++ ){
			r_[n][i] = 1.0f;
			pos_[n][i] = q_ / static_cast<float>( LW1 );		// p(z_n = i), i > 0
		}
		pos_[n][0] = 1 - q_;

		// when p(z_n > 0)
		// ij = i+j runs over all positions in sequence
		for( int ij = 0; ij < L; ij++ ){

			// extract (k+1)-mer y from positions (i-k,...,i)
			int y = posSeqs_[n]->extractKmer( ij,( ij < K ) ? ij : K );

			// j runs over all motif positions
			for( int j = ( 0 > ( ij-L+W ) ? 0 : ij-L+W ); j < ( W < (ij+1) ? W : ij+1 ); j++ ){

				// skip 'N' and other unknown alphabets
				if( y != -1 ){
					r_[n][L-W+1-ij+j] *= s_[y][j];
				} else {
					r_[n][L-W+1-ij+j] = 0.0f;
					break;
				}
			}
		}

		// calculate complete responsibilities and sum them up
		for( int i = 0; i < LW2; i++ ){
			if( r_[n][i] != 0.0f ){
				r_[n][i] *= pos_[n][i];
			}
			normFactor += r_[n][i];
		}

		// normalize responsibilities
		for( int i = 0; i < LW2; i++ ){
			r_[n][i] /= normFactor;
		}
		llikelihood_ += logf( normFactor );
	}
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
		int LW1 = L - W + 1;

		// ij = i+j runs over all positions in x
		for( int ij = 0; ij < L; ij++ ){
			int y = posSeqs_[n]->extractKmer( ij, ( ij < K ) ? ij : K );
			for( int j = ( 0 > ( ij-L+W ) ? 0 :  ij-L+W ); j < ( W < (ij+1) ? W : ij+1 ); j++ ){
				// skip 'N' and other unknown alphabets
				if( y != -1 && ( ij-j ) < LW1 ){
					n_[K][y][j] += r_[n][L-W+1-ij+j];
				}
			}
		}
	}

	// compute fractional occurrence counts from higher to lower order
	// k runs over all orders
	for( int k = K; k > 0; k-- ){
		for( int y = 0; y < Y_[k+1]; y++ ){
			int y2 = y % Y_[k];									// cut off the first nucleotide in (k+1)-mer
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
	bool iterate = true;							// flag for iterating before convergence
	int iteration = 0;

	int K = Global::modelOrder;
	int W = motif_->getW();

	int k, y, y2, j, i;

	float eta = Global::eta;					// learning rate for alpha learning

	// ToDo: write down log posterior and model v after each alpha updating step
	std::string opath = std::string( Global::outputDirectory ) + '/'
						+ Global::posSequenceBasename;
	std::string opath_log = opath + ".logposterior";
	std::string opath_vdiff = opath + ".vdiff";
	std::ofstream ofile_log( opath_log.c_str() );
	std::ofstream ofile_vdiff( opath_vdiff.c_str() );

	float lposterior_prev = 0.0f;
	float lposterior_new = 0.0f;
	float lposterior_diff = 0.0f;

	float*** v_prev = ( float*** )calloc( K+1, sizeof(float**) );
	for( k = 0; k < K+1; k++ ){
		v_prev[k] = ( float** )calloc( Y_[k+1], sizeof( float* ) );
		for( y = 0; y < Y_[k+1]; y++ ){
			v_prev[k][y] = ( float* )calloc( W, sizeof( float ) );
		}
	}

	// count the k-mers
	// 1. reset n_z_[k][y][j]
	for( k = 0; k < K+1; k++ ){
		for( y = 0; y < Y_[k+1]; y++ ){
			for( j = 0; j < W; j++ ){
				n_z_[k][y][j] = 0;
			}
		}
	}
	// 2. count k-mers for the highest order K
	// note: the last sequence is not counted, since during the first step
	// of iteration, it is added back.
	for( i = 0; i < ( int )posSeqs_.size()-1; i++ ){
		for( j = 0; j < W; j++ ){
			if( z_[i] > 0 ){
				y = posSeqs_[i]->extractKmer( z_[i]-1+j, ( z_[i]-1+j < K ) ?  z_[i]-1+j : K );
				if( y >= 0 ){
					n_z_[K][y][j]++;
				}
			}
		}
	}
	// 3. calculate nz for lower order k
	for( k = K; k > 0; k-- ){						// k runs over all lower orders
		for( y = 0; y < Y_[k+1]; y++ ){
			y2 = y % Y_[k];							// cut off the first nucleotide in (k+1)-mer
			for( j = 0; j < W; j++ ){
				n_z_[k-1][y2][j] += n_z_[k][y][j];
			}
		}
	}

	// write down the two terms of eq. 49 for plotting
	std::string opath_first = opath + "_debug.first";
	std::string opath_second = opath + "_debug.second";
	std::ofstream ofile_first( opath_first );
	std::ofstream ofile_second( opath_second );

	// iterate over
	while( iterate && iteration < Global::maxCGSIterations ){

		iteration++;

		std::cout << std::setprecision( 6 );
		std::cout << iteration << " iter:	" << std::endl;
		// get the previous v before updating alphas and z and q
		for( k = 0; k < K+1; k++ ){
			for( y = 0; y < Y_[k+1]; y++ ){
				for( j = 0; j < W; j++ ){
					v_prev[k][y][j] = motif_->getV()[k][y][j];
				}
			}
		}

		// sampling z and q
//		if( !Global::noZQSampling )		Gibbs_sampling_z_q();

		// update alphas:
		if( !Global::noAlphaUpdating ){

//			if( iteration % Global::interval == 0 ){

				lposterior_prev = calc_logPosterior_alphas( alpha_, K );

				// update alphas by stochastic optimization
				stochastic_optimize_alphas( K, W, eta, iteration );

				// update model parameter v
				motif_->updateVz_n( n_z_, alpha_, K );

				// todo: write out the difference of v's for comparison
				for( k = 0; k < K+1; k++ ){
					float v_diff = 0.0f;
					for( y = 0; y < Y_[k+1]; y++ ){
						for( j = 0; j < W; j++ ){
							v_diff += fabsf( motif_->getV()[k][y][j] - v_prev[k][y][j] );
						}
					}
					ofile_vdiff << std::setprecision(6) << v_diff << '\t';
				}
				ofile_vdiff << std::endl;
/*
				// update alphas by fixed learning rate
				for( int k = 0; k < K+1; k++ ){
					for( int j = k; j < W; j++ ){
						// update alphas:
						alpha_[k][j] = alpha_[k][j] + eta[k][j] * calc_gradient_alphas( alpha_, k, j );

						// update learning rates:
						eta[k][j] *= 0.9f;
					}
				}
*/

				lposterior_new = calc_logPosterior_alphas( alpha_, K );
//				lposterior_diff = lposterior_new - lposterior_prev;

//				std::cout << "old lpos= " << lposterior_prev << '\t';
//				std::cout << "new lpos= " << lposterior_new << '\t';
//				std::cout << "diff= " << lposterior_diff << std::endl;

				// write the log posterior into a file
				ofile_log << lposterior_new << std::endl;

//			}

		}

		if( Global::debugAlphas ){

			// update alphas
			for( k = 0; k < K+1; k++ ){
				for( j = 0; j < W; j++ ){
					alpha_[k][j] = ( float )iteration;
				}
			}

			// update v's
			motif_->updateVz_n( n_z_, alpha_, K );

			std::vector<float> first_term = debug_first_term_of_derivative( K );
			std::vector<float> second_term = debug_second_term_plus_prior_of_derivative( K );

			for( j = 0; j < W; j++ ){
				ofile_first << first_term[j] << '\t';
				ofile_second << second_term[j] << '\t';
			}
			ofile_first << std::endl;
			ofile_second << std::endl;
		}


		// only for writing out model after each iteration:
/*		motif_->calculateP();
		motif_->write( CGSIterations_ );*/
	}

	// obtaining a motif model

	// calculate probabilities
	motif_->calculateP();

	fprintf( stdout, "\n--- Runtime for Collapsed Gibbs sampling: %.4f seconds ---\n", ( ( float )( clock() - t0 ) ) / CLOCKS_PER_SEC );

	// free the memory:
	for( k = 0; k < K+1; k++ ){
		for( y = 0; y < Y_[k+1]; y++ ){
			free( v_prev[k][y] );
		}
		free( v_prev[k] );
	}
	free( v_prev );

}

void ModelLearning::Gibbs_sampling_z_q(){

	int N = ( int )posSeqs_.size();
	int W = motif_->getW();
	int K = Global::modelOrder;
	int K_bg = ( Global::bgModelOrder < K ) ? Global::bgModelOrder : K;
	int k, y, y_prev, y_bg, j, i, LW1, n, n_prev;
	int N_0 = 0;								// counts of sequences which do not contain motifs.

	float** v_bg = bg_->getV();

	// generated a random number, which is seeded by 42
	std::mt19937 rng( 42 );

	// sampling z:
	// loop over all sequences and drop one sequence each time and update r
	for( n = 0; n < N; n++ ){

		LW1 = posSeqs_[n]->getL() - W + 1;

		// get the index of the previous sequence
		n_prev = ( n > 0) ? n-1 : N-1;

		// calculate positional prior:
		pos_[n][0] = 1.0f - q_;
		for( i = 1; i <= LW1; i++ ){
			pos_[n][i] = q_ / ( float )LW1;
		}

		// count k-mers at position z[i]+j except the n'th sequence
		for( k = 0; k < K+1; k++ ){
			for( j = 0; j < W; j++ ){
				if( z_[n] != 0 ){
					// remove the k-mer counts from the current sequence with old z
					y = posSeqs_[n]->extractKmer( z_[n]-1+j, ( z_[n]-1+j < k ) ? z_[n]-1+j : k );
					if( y >= 0 ){
						n_z_[k][y][j]--;
					}
				}

				if( z_[n_prev] != 0 ){
					// add the k-mer counts from the previous sequence with updated z
					y_prev = posSeqs_[n_prev]->extractKmer( z_[n_prev]-1+j, ( z_[n_prev]-1+j < k ) ? z_[n_prev]-1+j : k );
					if( y_prev >= 0 ){
						n_z_[k][y_prev][j]++;
					}
				}
			}
		}

		// updated model parameters v excluding the n'th sequence
		motif_->updateVz_n( n_z_, alpha_, K );

		// compute log odd scores s[y][j], log likelihoods of the highest order K
		for( y = 0; y < Y_[K+1]; y++ ){
			for( j = 0; j < W; j++ ){
				y_bg = y % Y_[K_bg+1];
				s_[y][j] = motif_->getV()[K][y][j] / v_bg[K_bg][y_bg];
			}
		}

		// sampling equation:
		// calculate responsibilities over all LW1 positions on n'th sequence
		std::vector<float> posteriors;
		float normFactor = 0.0f;
		for( i = 1; i <= LW1; i++ ){
			r_[n][i] = 1.0f;
			for( j = 0; j < W; j++ ){
				// extract k-mers on the motif at position i over W of the n'th sequence
				y = posSeqs_[n]->extractKmer( i-1+j, ( i-1+j < K ) ? i-1+j : K );
				if( y >= 0 ){
					r_[n][i] *= s_[y][j];
				}
			}
			r_[n][i] *= pos_[n][i];
			normFactor += r_[n][i];
		}
		// for sequences that do not contain motif
		r_[n][0] = pos_[n][0];
		normFactor += r_[n][0];
		for( i = 0; i <= LW1; i++ ){
			r_[n][i] /= normFactor;
			posteriors.push_back( r_[n][i] );
		}
		// draw a new position z from discrete posterior distribution
		std::discrete_distribution<> posterior_dist( posteriors.begin(), posteriors.end() );
		z_[n] = posterior_dist( rng );					// draw a sample z randomly
		if( z_[n] == 0 ) N_0++;
	}

	// sampling q:
	// draw two random numbers Q and P from Gamma distribution
	std::gamma_distribution<> P_Gamma_dist( N_0 + 1, 1 );
	std::gamma_distribution<> Q_Gamma_dist( N - N_0 + 1, 1 );
	double P = P_Gamma_dist( rng );						// draw a sample for P
	double Q = Q_Gamma_dist( rng );						// draw a sample for Q

	q_ = ( float ) Q / ( float )( Q + P );				// calculate q_

}

void ModelLearning::stochastic_optimize_alphas( int K, int W, float eta, int iter ){
	// update alphas using stochastic optimization algorithm ADAM (DP Kingma & JL Ba 2015)

	int k, j;
	float beta1 = 0.9f;									// exponential decay rate for the moment estimates
//	float beta2 = 0.999f;								// exponential decay rate for the moment estimates
	float beta2 = 0.99f;
	float epsilon = 1e-8f;								// cutoff
	float gradient;										// gradient of log posterior of alpha
	float m1;											// first moment vector (the mean)
	float m2;											// second moment vector (the uncentered variance)

	float t = ( float )iter;

	for( k = 0; k < K+1; k++ ){

		for( j = 0; j < W; j++ ){
			// reparameterise alpha on log scale: alpha = e^a
			float a  = logf( alpha_[k][j] );

			// get gradients w.r.t. stochastic objective at timestep t
			gradient = alpha_[k][j] * calc_gradient_alphas( alpha_, k, j );

			// update biased first moment estimate
			m1_t_[k][j] = beta1 * m1_t_[k][j] + ( 1 - beta1 ) * gradient;

			// update biased second raw moment estimate
			m2_t_[k][j] = beta2 * m2_t_[k][j] + ( 1 - beta2 ) * gradient * gradient;

			// compute bias-corrected first moment estimate
			m1 = m1_t_[k][j] / ( 1 - powf( beta1, t ) );

			// compute bias-corrected second raw moment estimate
			m2 = m2_t_[k][j] / ( 1 - powf( beta2, t ) );

			// update parameter a due to alphas
			a -= eta * m1 / ( ( sqrtf( m2 ) + epsilon ) * sqrtf( t ) );

			alpha_[k][j] = expf( a );
		}

	}

}

float ModelLearning::calc_logPosterior_alphas( float** alpha, int k ){
	// calculate partial log posterior of alphas due to equation 46 in the theory

	float logPosterior = 0.0f;
	float*** v = motif_->getV();
	float** v_bg = bg_->getV();
	int W = motif_->getW();
	int N = static_cast<int>( posSeqs_.size() ) - 1;
	int j, y, y2;

	for( j = 0; j < W; j++ ){

		// the first term of equation 46
		logPosterior -= 2.0f * logf( alpha[k][j] );

		// the second term of equation 46
		logPosterior -= Global::modelBeta * powf( Global::modelGamma, ( float )k ) / alpha[k][j];

		// the third term of equation 46
		logPosterior += ( float )ipow( Y_[1], k ) * boost::math::lgamma( alpha[k][j] );

		// the forth term of equation 46
		for( y = 0; y < Y_[k+1]; y++ ){
			y2 = y % Y_[k];									// cutoff the first nucleotide in the (k+1)-mer
			if( k == 0 ){
				// the first part
				logPosterior += boost::math::lgamma( ( float )n_z_[k][y][j] + alpha[k][j] * v_bg[k][y] );

				// the second part
				logPosterior -= boost::math::lgamma( alpha[k][j] * v_bg[k][y] );

			} else {
				// the first part
				logPosterior += boost::math::lgamma( ( float )n_z_[k][y][j] + alpha[k][j] * v[k-1][y2][j] );

				// the second part
				logPosterior -= boost::math::lgamma( alpha[k][j] * v[k-1][y2][j] );

			}
		}

		// the last term of equation 46
		for( y = 0; y < Y_[k]; y++ ){
			if( j == 0 ){
				logPosterior -= boost::math::lgamma( ( float )N + alpha[k][j] );
			} else {
				logPosterior -= boost::math::lgamma( ( float )n_z_[k][y][j-1] + alpha[k][j] );
			}
		}
	}

	return logPosterior;
}

float ModelLearning::calc_gradient_alphas( float** alpha, int k, int j ){
	// calculate partial gradient of the log posterior of alphas due to equation 47 in the theory

	float gradient = 0.0f;
	float*** v = motif_->getV();
	float** v_bg = bg_->getV();
	int N = static_cast<int>( posSeqs_.size() ) - 1;
	int y, y2;


	// the first term of equation 47, part of the prior
//	gradient -= 2.0f / alpha[k][j];

	// the second term of equation 47, part of the prior
//	gradient += Global::modelBeta * powf( Global::modelGamma, ( float )k ) / powf( alpha[k][j], 2.f );

	// the third term of equation 47
	gradient += ( float )ipow( Y_[1], k ) * boost::math::digamma( alpha[k][j] );


	// the forth term of equation 47
	for( y = 0; y < Y_[k+1]; y++ ){
		y2 = y % Y_[k];
		if( k == 0 ){
			// the first part
			gradient += v_bg[k][y] * boost::math::digamma( ( float )n_z_[k][y][j] + alpha[k][j] * v_bg[k][y] );

			// the second part
			gradient -= v_bg[k][y] * boost::math::digamma( alpha[k][j] * v_bg[k][y] );

		} else {
			// the first part
			gradient += v[k-1][y2][j] * boost::math::digamma( ( float )n_z_[k][y][j] + alpha[k][j] * v[k-1][y2][j] );

			// the second part
			gradient -= v[k-1][y2][j] * boost::math::digamma( alpha[k][j] * v[k-1][y2][j] );

		}
	}

	// the last term of equation 47
	for( y = 0; y < Y_[k]; y++ ){

		if( j == 0 ){

			gradient -= boost::math::digamma( ( float )N + alpha[k][j] );

		} else {

			gradient -= boost::math::digamma( ( float )n_z_[k][y][j-1] + alpha[k][j] );
		}
	}

	return gradient;
}

void ModelLearning::debug_optimization_alphas( float** alphas, int K, int W ){
	// test if the log posterior function of alpha fits to its gradient

	int k, j;

	for( k = 0; k < K+1; k++ ){

		std::cout << "k = " << k << std::endl;

		float stepsize = 0.0001f;
		float diff;
		for( j = 0; j < W; j++ ){
			std::cout << "j = " << j << '\t';
			std::cout << std::setprecision(8) << "alpha = " << alphas[k][j] << '\t';
			// calculate the delta of log posterior
			alphas[k][j] += stepsize;
			diff = calc_logPosterior_alphas( alphas, k );

			alphas[k][j] -= 2.0f * stepsize;
			diff -= calc_logPosterior_alphas( alphas, k );

			diff /= ( 2.0f * stepsize );

			// recover alphas
			alphas[k][j] += stepsize;

			// print out the gradient

			std::cout << "analytical gradient: "<< calc_gradient_alphas( alphas, k, j ) << std::endl;
	//		std::cout << "numerical gradient: " << std::setprecision(12) << diff << std::endl;

		}

	}
	std::cout << std::endl;
}

std::vector<float> ModelLearning::debug_first_term_of_derivative( int k ){

	float*** v = motif_->getV();
	std::vector<float> sum;

	size_t count = 0;
	for( int j = 0; j < motif_->getW(); j++ ){
		float s = 0.0f;
		for( int y = 0; y < Y_[k+1]; y++ ){
			count++;
			int y2 = y % Y_[k];
			std::cout << count << '\t'
					<< v[k-1][j][y2] << '\t'
					<< v[k][y][j] << '\t'
					<< v[k-1][y2][j] << '\t';
			s += v[k-1][j][y2] * ( logf( v[k][y][j] / v[k-1][y2][j] ) );
			std::cout << s << std::endl;
		}
		sum.push_back( s );
	}

	return sum;

}

std::vector<float> ModelLearning::debug_second_term_plus_prior_of_derivative( int k ){

	float*** v = motif_->getV();
	std::vector<float> sum;

	float N = ( float )posSeqs_.size();

	for( int j = 0; j < motif_->getW(); j++ ){
		float s = 0.0f;
		// first sum up prior
		float prior = -2.0f / alpha_[k][j] + Global::modelBeta *
				powf( Global::modelGamma, ( float )k ) / powf( alpha_[k][j], 2.0f );
		float first_term = 0.0f;
		float second_term = 0.0f;
		for( int y = 0; y < Y_[k]; y++ ){
			if( j > 0 ){
				first_term = 1.5f * ( float )n_z_[k][y][j-1] / ( ( ( float )n_z_[k][y][j-1] + alpha_[k][j] ) * alpha_[k][j] );
			} else {
				first_term = 1.5f * ( float )n_z_[k][y][j-1] / ( ( N + alpha_[k][j] ) * alpha_[k][j] );
			}
		}
		for( int y = 0; y < Y_[k+1]; y++ ){
			int y2 = y % Y_[k];
			int yk = y / Y_[1];
			if( j > 0 ){
				second_term = 0.5f * ( v[k][y][j] - v[k-1][y2][j] ) / ( v[k][y][j] * ( ( float )n_z_[k][yk][j-1] + 1.0e-8f ) );
			} else {
				second_term = 0.5f * ( v[k][y][j] - v[k-1][y2][j] ) / ( v[k][y][j] * N );
			}
		}
		s = prior + first_term + second_term;

		sum.push_back( s );
	}
	return sum;
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

	if( Global::EM ){
		// output (k+1)-mer counts n[k][y][j]
		std::string opath_n = opath + ".EMcount";
		std::ofstream ofile_n( opath_n.c_str() );
		for( j = 0; j < W; j++ ){
			for( k = 0; k < K+1; k++ ){
				for( y = 0; y < Y_[k+1]; y++ ){
					ofile_n << ( int )n_[k][y][j] << ' ';
				}
				ofile_n << std::endl;
			}
			ofile_n << std::endl;
		}

		// output responsibilities r[n][i] and position(s) of motif(s) pos_[n][i]
		std::string opath_r = opath + ".EMweight";
		std::string opath_pos = opath + ".EMposition";
		std::ofstream ofile_r( opath_r.c_str() );
		std::ofstream ofile_pos( opath_pos.c_str() );
		ofile_pos << "seq" << '\t' << "positions" << std::endl;
		float cutoff = 0.3f;									// threshold for having a motif on the sequence in term of responsibilities
		for( size_t n = 0; n < posSeqs_.size(); n++ ){
			int LW1 = posSeqs_[n]->getL() - W + 1;
			ofile_r << std::scientific << std::setprecision( 2 ) << r_[n][0] << ' ';	// print out the responsibility of not having a motif on the sequence
			ofile_pos << n+1 << '\t';							// print out the sequence number
			for( i = LW1; i > 0; i-- ){
				ofile_r << std::setprecision( 2 ) << r_[n][i] << ' ';
				if( r_[n][i] >= cutoff ){
					ofile_pos << LW1-i+1 << '\t';
				}
			}
			ofile_r << std::endl;
			ofile_pos << std::endl;
		}

		// output parameter alphas alpha[k][j]
		std::string opath_alpha = opath + ".EMalpha";
		std::ofstream ofile_alpha( opath_alpha.c_str() );
		for( k = 0; k < K+1; k++ ){
			ofile_alpha << "k = " << k << std::endl;
			for( j = 0; j < W; j++ ){
				ofile_alpha << std::setprecision( 3 ) << alpha_[k][j] << ' ';
			}
			ofile_alpha << std::endl;
		}

	} else if( Global::CGS ){
		// output (k+1)-mer integral counts nz[k][y][j]
		std::string opath_n = opath + ".CGScount";
		std::ofstream ofile_n( opath_n.c_str() );
		for( j = 0; j < W; j++ ){
			for( k = 0; k < K+1; k++ ){
				for( y = 0; y < Y_[k+1]; y++ ){
					ofile_n << n_z_[k][y][j] << ' ';
				}
				ofile_n << std::endl;
			}
			ofile_n << std::endl;
		}

		// output responsibilities r[n][i]
		std::string opath_r = opath + ".CGSweight";
		std::ofstream ofile_r( opath_r.c_str() );
		for( size_t n = 0; n < posSeqs_.size(); n++ ){
			for( i = 0; i < posSeqs_[n]->getL()-W+2; i++ ){
				ofile_r << std::scientific << std::setprecision( 2 ) << r_[n][i] << ' ';
			}
			ofile_r << std::endl;
		}

		// output parameter alphas alpha[k][j]
		std::string opath_alpha = opath + ".CGSalpha";
		std::ofstream ofile_alpha( opath_alpha.c_str() );
		for( k = 0; k < K+1; k++ ){
			ofile_alpha << "k = " << k << std::endl;
			for( j = 0; j < W; j++ ){
				ofile_alpha << std::setprecision( 3 ) << alpha_[k][j] << ' ';
			}
			ofile_alpha << std::endl;
		}

		// output positions of motifs z_[n]
		std::string opath_z = opath + ".CGSposition";
		std::ofstream ofile_z( opath_z.c_str() );
		ofile_z << "seq" << '\t' << "start" <<'\t' << "end" << std::endl;
		for( size_t n = 0; n < posSeqs_.size(); n++ ){
			ofile_z << n+1 << '\t' << z_[n] <<'\t' << z_[n]+W-1 << std::endl;
		}
	}
}
