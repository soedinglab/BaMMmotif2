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

ModelLearning::ModelLearning( Motif* motif, BackgroundModel* bg, std::vector<size_t> folds ){

	motif_ = motif;
	bg_ = bg;
	K_ = motif_->getK();
	W_ = motif_->getW();
	q_ = Global::q;

	// define parameters
	for( size_t k = 0; k < K_+8; k++ ){
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
	z_ = ( size_t* )calloc( posSeqs_.size(), sizeof( int ) );

	for( size_t n = 0; n < posSeqs_.size(); n++ ){
		size_t LW2 = posSeqs_[n]->getL() - W_ + 2;
		r_[n] = ( float* )calloc( LW2, sizeof( float ) );
		pos_[n] = ( float* )calloc( LW2, sizeof( float ) );
	}

	// allocate memory for n_[k][y][j] and probs_[k][y][j] and initialize them
	n_ = ( float*** )calloc( K_+1, sizeof( float** ) );
	n_z_ = ( size_t*** )calloc( K_+1, sizeof( size_t** ) );
	for( size_t k = 0; k < K_+1; k++ ){
		// allocate -1 position for y
		n_[k] = ( float** )calloc( Y_[k+1]+1, sizeof( float* ) )+1;
		n_z_[k] = ( size_t** )calloc( Y_[k+1]+1, sizeof( size_t* ) )+1;
		// allocate -K positions for j
		for( int y = -1; y < static_cast<int>( Y_[k+1] ); y++ ){
			n_[k][y] = ( float* )calloc( W_+K_, sizeof( float ) )+K_;
			n_z_[k][y] = ( size_t* )calloc( W_+K_, sizeof( size_t ) )+K_;
		}
	}

	// allocate memory for alpha_[k][j] and initialize it
	alpha_ = ( double** )malloc( ( K_+1 ) * sizeof( double* ) );
	m1_t_ = ( double** )calloc( K_+1, sizeof( double* ) );
	m2_t_ = ( double** )calloc( K_+1, sizeof( double* ) );
	for( size_t k = 0; k < K_+1; k++ ){
		alpha_[k] = ( double* )malloc( W_ * sizeof( double ) );
		m1_t_[k] = ( double* )calloc( W_, sizeof( double ) );
		m2_t_[k] = ( double* )calloc( W_, sizeof( double ) );
		for( size_t j = 0; j < W_; j++ ){
			alpha_[k][j] = Global::modelAlpha[k];
		}
	}

}

ModelLearning::~ModelLearning(){

	for( size_t n = 0; n < posSeqs_.size(); n++ ){
		free( r_[n] );
		free( pos_[n] );
	}
	free( r_ );
	free( pos_ );
	free( z_ );

	for( size_t k = 0; k < K_+1; k++ ){
		for( int y = -1; y < static_cast<int>( Y_[k+1] ); y++ ){
			free( n_[k][y] - K_ );
			free( n_z_[k][y] - K_ );
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

	if( Global::verbose ){
		fprintf( stderr," ______\n"
						"|      |\n"
						"|  EM  |\n"
						"|______|\n\n" );
	}

	clock_t t0 = clock();
	bool 	iterate = true;									// flag for iterating before convergence

	float 	v_diff, llikelihood_prev, llikelihood_diff = 0.0f;
	float**	v_before;										// hold the parameters of the highest-order before EM

	// allocate memory for parameters v[y][j] with the highest order
	v_before = ( float** )calloc( Y_[K_+1], sizeof( float* ) );
	for( size_t y = 0; y < Y_[K_+1]; y++ ){
		v_before[y] = ( float* )calloc( W_, sizeof( float ) );
	}

	int EMIterations = 0;
	// iterate over
	while( iterate && ( EMIterations < Global::maxEMIterations ) ){

		EMIterations++;

		// get parameter variables with highest order before EM
		llikelihood_prev = llikelihood_;
		for( size_t y = 0; y < Y_[K_+1]; y++ ){
			for( size_t j = 0; j < W_; j++ ){
				v_before[y][j] = motif_->getV()[K_][y][j];
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
		for( size_t y = 0; y < Y_[K_+1]; y++ ){
			for( size_t j = 0; j < W_; j++ ){
				v_diff += fabsf( motif_->getV()[K_][y][j] - v_before[y][j] );
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
	for( size_t y = 0; y < Y_[K_+1]; y++ ){
		free( v_before[y] );
	}
	free( v_before );

	fprintf( stdout, "\n--- Runtime for EM: %.4f seconds ---\n", ( ( float )( clock() - t0 ) ) / CLOCKS_PER_SEC );

    return 0;
}

void ModelLearning::EStep(){

	llikelihood_ = 0.0f;

	motif_->calculateLinearS( bg_->getV() );

	float** s = motif_->getS();
//	int N0 = 0;
	// parallel the code
//	#pragma omp parallel for

	// calculate responsibilities r_[n][i] at position i in sequence n
	// n runs over all sequences
	for( size_t n = 0; n < posSeqs_.size(); n++ ){

		size_t L = posSeqs_[n]->getL();
		size_t LW1 = L - W_ + 1;
		size_t LW2 = L - W_ + 2;
		int* kmer = posSeqs_[n]->getKmer();
		float normFactor = 0.0f;

		// reset r_[n][i] and pos_[n][i]
		for( size_t i = 0; i < LW2; i++ ){
			r_[n][i] = 1.0f;
			pos_[n][i] = q_ / static_cast<float>( LW1 );
		}
		pos_[n][0] = 1 - q_;

		// when p(z_n > 0)
		// ij = i+j runs over all positions in sequence
		for( size_t ij = 0; ij < L; ij++ ){

			// extract (K+1)-mer y from positions (i-k,...,i)
			int y = ( kmer[ij] >= 0 ) ? kmer[ij] % static_cast<int>( Y_[K_+1] ) : -1;
			// j runs over all motif positions
			size_t padding = ( static_cast<int>( ij-L+W_ ) > 0 ) * ( ij-L+W_ );
			for( size_t j = padding; j < ( W_ < (ij+1) ? W_ : ij+1 ); j++ ){
				r_[n][LW1-ij+j] *= s[y][j];
			}
		}

		// calculate complete responsibilities and sum them up
		for( size_t i = 0; i < LW2; i++ ){
			r_[n][i] *= pos_[n][i];
			normFactor += r_[n][i];
		}

		// normalize responsibilities
		for( size_t i = 0; i < LW2; i++ ){
			r_[n][i] /= normFactor;
		}

//		if( r_[n][0] > 0.7 ) N0++;
		// calculate log likelihood over all sequences
		llikelihood_ += logf( normFactor );
	}

//	std::cout << "N0=" << N0 << std::endl;
}

void ModelLearning::MStep(){

	// reset the fractional counts n
	for( size_t k = 0; k < K_+1; k++ ){
		for( size_t y = 0; y < Y_[k+1]; y++ ){
			for( size_t j = 0; j < W_; j++ ){
				n_[k][y][j] = 0.0f;
			}
		}
	}

	// compute fractional occurrence counts for the highest order K
	// n runs over all sequences
	for( size_t n = 0; n < posSeqs_.size(); n++ ){
		size_t L = posSeqs_[n]->getL();
		int* kmer = posSeqs_[n]->getKmer();

		// ij = i+j runs over all positions i on sequence n
		for( size_t ij = 0; ij < L; ij++ ){

			int y = ( kmer[ij] >= 0 ) ? kmer[ij] % static_cast<int>( Y_[K_+1] ) : -1;

			size_t padding = ( static_cast<int>( ij-L+W_ ) > 0 ) * ( ij-L+W_ );

			for( size_t j = padding; j < ( W_ < (ij+1) ? W_ : ij+1 ); j++ ){

				n_[K_][y][j] += r_[n][L-W_+1-ij+j];

			}
		}
	}

	// compute fractional occurrence counts from higher to lower order
	// k runs over all orders
	for( size_t k = K_; k > 0; k-- ){

		for( size_t y = 0; y < Y_[k+1]; y++ ){

			size_t y2 = y % Y_[k];

			for( size_t j = 0; j < W_; j++ ){

				n_[k-1][y2][j] += n_[k][y][j];

			}
		}
	}

	// update model parameters v[k][y][j]
	motif_->updateV( n_, alpha_, K_ );
}

void ModelLearning::optimize_q(){

	// optimize hyper-parameter q
	// motif.updateV()
}

void ModelLearning::GibbsSampling(){

	if( Global::verbose ){
		fprintf(stderr, " ___________________________\n"
						"|                           |\n"
						"|  Collapsed Gibbs sampler  |\n"
						"|___________________________|\n\n" );
	}
	clock_t t0 = clock();
	size_t iteration = 0;

	float eta = Global::eta;						// learning rate for alpha learning

	// initialize z for all the sequences
	if( !Global::noInitialZ ){

		// E-step: calculate posterior
		EStep();

		// extract initial z from the indices of the largest responsibilities
		for( size_t n = 0; n < posSeqs_.size(); n++ ){
			size_t LW1 = posSeqs_[n]->getL() - W_ + 1;
			float maxR = r_[n][0];
			size_t maxIdx = 0;
			for( size_t i = 1; i <= LW1; i++ ){
				if( r_[n][i] > maxR ){
					maxR = r_[n][i];
					maxIdx = LW1+1-i;
				}
			}
			z_[n] = maxIdx;
		}

	} else {
		// initialize z with a random number
		for( size_t n = 0; n < posSeqs_.size(); n++ ){
			size_t LW2 = posSeqs_[n]->getL() - motif_->getW() + 2;
			z_[n] = ( size_t )rand() % LW2;
		}
	}

	// count the k-mers
	// 1. reset n_z_[k][y][j] to 0
	for( size_t k = 0; k < K_+1; k++ ){
		for( size_t y = 0; y < Y_[k+1]; y++ ){
			for( int j = -(int)K_; j < (int)W_; j++ ){
				n_z_[k][y][j] = 0;
			}
		}
	}

	// 2. count k-mers for the highest order K
	for( size_t n = 0; n < posSeqs_.size(); n++ ){
		if( z_[n] > 0 ){
			int* kmer = posSeqs_[n]->getKmer();
			for( int j = ( z_[n] <= K_ ) ? 1-(int)z_[n] : -(int)K_; j < (int)W_; j++ ){
				int pos = static_cast<int>( z_[n] )-1+j;
				int y = ( kmer[pos] >= 0 ) ? kmer[pos] % static_cast<int>( Y_[K_+1] ) : -1;
				n_z_[K_][y][j]++;
			}
		}
	}

	// compute k-mer counts for all the lower orders
	for( size_t k = K_; k > 0; k-- ){
		for( size_t y = 0; y < Y_[k+1]; y++ ){
			size_t y2 = y % Y_[k];
			for( int j = -(int)K_; j < (int)W_; j++ ){
				n_z_[k-1][y2][j] += n_z_[k][y][j];
			}
		}
	}

	// Gibbs sampling position z and fraction q
	for( size_t iter = 0; iter < 10; iter++ ){
		Gibbs_sample_z_q();
	}

	// vector to store the last a few alphas for sampling methods
	std::vector<std::vector<double>> alpha_avg( K_+1 );
	for( size_t k = 0; k < K_+1; k++ ){
		alpha_avg[k].resize( W_ );
		for( size_t j = 0; j < W_; j++ ){
			alpha_avg[k][j] = 0.0;
		}
	}

	// iterate over
	while( iteration < Global::maxCGSIterations ){

		iteration++;

		// Gibbs sampling position z and fraction q
		if( !Global::noZSampling )	Gibbs_sample_z_q();

		// update alphas by stochastic optimization
		if( !Global::noAlphaOptimization ){

			stochastic_optimize_alphas( K_, W_, eta, iteration );

		} else if( Global::GibbsMHalphas ){

			GibbsMH_sample_alphas( iteration );

			if( iteration > Global::maxCGSIterations - 10 ){
				for( size_t k = 0; k < K_+1; k++ ){
					for( size_t j = 0; j < W_; j++ ){
						alpha_avg[k][j] += alpha_[k][j];
					}
				}
			}

		} else if( Global::dissampleAlphas ){

			discrete_sample_alphas( iteration );

			if( iteration > Global::maxCGSIterations - 10 ){
				for( size_t k = 0; k < K_+1; k++ ){
					for( size_t j = 0; j < W_; j++ ){
						alpha_avg[k][j] += alpha_[k][j];
					}
				}
			}
		}
	}

	// obtaining a motif model:
	if( Global::GibbsMHalphas || Global::dissampleAlphas ){
		// average alphas over the last few steps for GibbsMH
		for( size_t k = 0; k < K_+1; k++ ){
			for( size_t j = 0; j < W_; j++ ){
				alpha_[k][j] = alpha_avg[k][j] / 10.0;
			}
		}
	}

	// update model parameter v
	motif_->updateVz_n( n_z_, alpha_, K_ );

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

	size_t N = posSeqs_.size();
	size_t K_bg = ( Global::bgModelOrder < K_ ) ? Global::bgModelOrder : K_;
	size_t N_0 = 0;
	float N_1 = static_cast<float>( posSeqs_.size() ) * q_;

	llikelihood_ = 0.0f;

	float*** v = motif_->getV();
	float** v_bg = bg_->getV();

	// updated model parameters v excluding the n'th sequence
	motif_->updateVz_n( n_z_, alpha_, K_ );
	// compute log odd scores s[y][j], log likelihoods of the highest order K
	motif_->calculateLinearS( v_bg );

	float** s = motif_->getS();

	// sampling z:
	// loop over all sequences and drop one sequence each time and update r
	for( size_t n = 0; n < N; n++ ){

		size_t L = posSeqs_[n]->getL();
		size_t LW1 = L - W_ + 1;
		int* kmer = posSeqs_[n]->getKmer();

		// count k-mers at position z_i+j except the n'th sequence
		// remove the k-mer counts from the sequence with the current z
		// and re-calculate the log odds scores in linear space
		/*
		 * -------------- faster version of removing k-mer ------------------
		 */

		if( z_[n] > 0 ){

			for( size_t k = 0; k < K_+1; k++ ){

				for( int j = ( z_[n] <= K_ ) ? 1-(int)z_[n] : -(int)K_; j < (int)W_; j++ ){
					int pos = static_cast<int>( z_[n] )-1+j;
					int y = ( kmer[pos] >= 0 ) ? kmer[pos] % static_cast<int>( Y_[k+1] ) : -1;

					n_z_[k][y][j]--;

					if( j >= 0 && y >= 0 ){
						if( k == 0 ){
							v[k][y][j]= ( static_cast<float>( n_z_[k][y][j] ) + static_cast<float>( alpha_[k][j] ) * v_bg[k][y] )
										/ ( N_1 + static_cast<float>( alpha_[k][j] ) );
						} else {

							size_t y2 = static_cast<size_t>( y ) % Y_[k];
							size_t yk = static_cast<size_t>( y ) / Y_[1];

							v[k][y][j] = ( static_cast<float>( n_z_[k][y][j] ) + static_cast<float>( alpha_[k][j] ) * v[k-1][y2][j] )
										/ ( static_cast<float>( n_z_[k-1][yk][j-1] ) + static_cast<float>( alpha_[k][j] ) );
						}

						if( k == K_ ){

							size_t y_bg = static_cast<size_t>( y ) % Y_[K_bg+1];

							s[y][j] = v[K_][y][j] / v_bg[K_bg][y_bg];

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
			for( j = ( z_[n] <= K ) ? 1-z_[n] : -K; j < W_; j++ ){
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
		for( size_t i = 1; i <= LW1; i++ ){
			pos_[n][i] = pos_i;
			r_[n][i] = 1.0f;
		}

		// todo: could be parallelized by extracting 8 sequences at once
		// ij = i+j runs over all positions in sequence
		for( size_t ij = 0; ij < L; ij++ ){

			// extract (K+1)-mer y from positions (i-k,...,i)
			int y = ( kmer[ij] >= 0 ) ? kmer[ij] % static_cast<int>( Y_[K_+1] ) : -1;

			// j runs over all motif positions
			size_t padding = ( static_cast<int>( ij-L+W_ ) > 0 ) * ( ij-L+W_ );
			for( size_t j = padding; j < ( W_ < (ij+1) ? W_ : ij+1 ); j++ ){
				r_[n][LW1-ij+j] *= s[y][j];
			}
		}
		for( size_t i = 1; i <= LW1; i++ ){
			r_[n][i] *= pos_[n][i];
			normFactor += r_[n][i];
		}

		// calculate log likelihood of sequences given the motif positions z and model parameter v
		llikelihood_ += logf( normFactor );

		// normalize responsibilities and append them to an array
		std::vector<float> posteriors;
		r_[n][0] /= normFactor;
		posteriors.push_back( r_[n][0] );

		for( size_t i = LW1; i >= 1; i-- ){
			r_[n][i] /= normFactor;
			posteriors.push_back( r_[n][i] );
		}

		// draw a new position z from the discrete distribution of posterior
		std::discrete_distribution<size_t> posterior_dist( posteriors.begin(), posteriors.end() );
		z_[n] = posterior_dist( Global::rngx );

		if( z_[n] == 0 ){
			// count sequences which do not contain motifs.
			N_0++;

		} else {
			// add the k-mer counts from the current sequence with the updated z
			for( int j = ( z_[n] <= K_ ) ? 1-(int)z_[n] : -(int)K_; j < (int)W_; j++ ){
				int pos = static_cast<int>( z_[n] )-1+j;
				for( size_t k = 0; k < K_+1; k++ ){
					int y = ( kmer[pos] >= 0 ) ? kmer[pos] % static_cast<int>( Y_[k+1] ) : -1;
					n_z_[k][y][j]++;
				}
			}
		}

	}

/*
	// sampling q:
	if( !Global::noQSampling ){
		boost::math::beta_distribution<> q_beta_dist( N - N_0 + 1, N_0 + 1);
		q_ = ( float )quantile( q_beta_dist, ( float )rand() / ( float )RAND_MAX );
	}
*/

}

void ModelLearning::stochastic_optimize_alphas( size_t K, size_t W_, float eta, size_t iter ){
	// update alphas using stochastic optimization algorithm ADAM (DP Kingma & JL Ba 2015)

	double beta1 = 0.9;		// exponential decay rate for the moment estimates
	double beta2 = 0.99;	// exponential decay rate for the moment estimates
	double epsilon = 1e-8;	// cutoff
	double gradient;		// gradient of log posterior of alpha
	double m1;				// first moment vector (the mean)
	double m2;				// second moment vector (the uncentered variance)

	double t = static_cast<double>( iter );

	for( size_t k = 0; k < K+1; k++ ){

		for( size_t j = 0; j < W_; j++ ){

			// re-parameterise alpha on log scale: alpha = e^a
			double a = log( alpha_[k][j] );

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

void ModelLearning::GibbsMH_sample_alphas( size_t iter ){
	// Gibbs sampling alphas in exponential space with Metropolis-Hastings algorithm

	for( size_t k = 0; k < K_+1; k++ ){

		for( size_t j = 0; j < W_; j++ ){

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

void ModelLearning::discrete_sample_alphas( size_t iter ){
	// sample an alpha from the discrete distribution of its log posterior

	for( size_t k = 0; k < K_+1; k++ ){

		for( size_t j = 0; j < W_; j++ ){

			std::vector<double> condProb;

			double condProb_new;

			double base = calc_logCondProb_a( iter, 0.0, k, j );

			for( size_t it = 0; it < 100; it++ ){

				condProb_new = exp( calc_logCondProb_a( iter, static_cast<double>( it ) / 10.0, k, j ) - base );

				condProb.push_back( condProb_new );

			}

			std::discrete_distribution<> posterior_dist( condProb.begin(), condProb.end() );

			alpha_[k][j] = exp( posterior_dist( Global::rngx ) / 10.0 );
		}
	}
}

double ModelLearning::calc_gradient_alphas( double** alpha, size_t k, size_t j ){
	// calculate partial gradient of the log posterior of alphas due to equation 47 in the theory
	// Note that j >= k

	double gradient = 0.0;
	float*** v = motif_->getV();
	float** v_bg = bg_->getV();
	size_t N = posSeqs_.size() - 1;

	// the first term of equation 47
	gradient -= 2.0 / alpha[k][j];

	// the second term of equation 47
	gradient += Global::modelBeta * pow( Global::modelGamma, static_cast<double>( k ) ) / pow( alpha[k][j], 2.0 );

	// the third term of equation 47
	gradient += static_cast<double>( ipow( Y_[1], k ) ) * boost::math::digamma( alpha[k][j] );

	// the forth term of equation 47
	for( size_t y = 0; y < Y_[k+1]; y++ ){

		if( k == 0 ){
			// the first part
			gradient += v_bg[k][y] * boost::math::digamma( static_cast<double>( n_z_[k][y][j] ) + alpha[k][j] * v_bg[k][y] );

			// the second part
			gradient -= v_bg[k][y] * boost::math::digamma( alpha[k][j] * v_bg[k][y] );

		} else {

			size_t y2 = y % Y_[k];

			// the first part
			gradient += v[k-1][y2][j] * boost::math::digamma( static_cast<double>( n_z_[k][y][j] ) + alpha[k][j] * v[k-1][y2][j] );

			// the second part
			gradient -= v[k-1][y2][j] * boost::math::digamma( alpha[k][j] * v[k-1][y2][j] );

		}

	}

	// the last term of equation 47
	for( size_t y = 0; y < Y_[k]; y++ ){

		if( k == 0 ){

			gradient -= boost::math::digamma( static_cast<double>( N ) + alpha[k][j] );

		} else if( j == 0 ){
			size_t sum = 0;
			for( size_t  a = 0; a < Y_[1]; a++ ){
				size_t  ya = y * Y_[1] + a;
				sum += n_z_[k][ya][j];
			}
			gradient -= boost::math::digamma( static_cast<double>( sum ) + alpha[k][j] );

		} else {

			gradient -= boost::math::digamma( static_cast<double>( n_z_[k-1][y][j-1] ) + alpha[k][j] );

		}
	}
	return gradient;
}

double ModelLearning::calc_logCondProb_a( size_t iteration, double a, size_t k, size_t j ){
	// calculate partial log conditional probabilities of a's due to equation 50 in the theory

	double logCondProbA = 0.0;
	float*** v = motif_->getV();
	float** v_bg = bg_->getV();
	size_t N = posSeqs_.size() - 1;
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

		for( size_t y = 0; y < Y_[k]; y++ ){

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
		for( size_t y = 0; y < Y_[k]; y++ ){

			// the third term of equation 50
			logCondProbA += boost::math::lgamma( alpha );

			if( be_printed ) std::cout << boost::math::lgamma( alpha ) << '\t' << logCondProbA << "\t(3+" << y <<") : lgamma(" << alpha << ")"<< std::endl;


			// the forth term of equation 50
			for( size_t A = 0; A < Y_[1]; A++ ){

				size_t ya = y * Y_[1] + A;

				size_t y2 = ya % Y_[k];

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
			if( n_z_[k-1][y][j-1] <= 1 ){

				logCondProbA = -a - Global::modelBeta * pow( Global::modelGamma, static_cast<double>( k ) ) / alpha;

				if( be_printed ) std::cout << "* 6 : " << logCondProbA << " since n_[k-1][y][j-1]<=1" << std::endl;

			}

			if( be_printed ) std::cout << std::endl;
		}

	}

	return logCondProbA;
}

double ModelLearning::calc_prior_alphas( double** alpha, size_t k ){
	// calculate partial log conditional probabilities of alphas due to equation 34 in the theory

	double logPrior = 0.0;
	for( size_t j = 0; j < W_; j++ ){

		// the first term of equation 46
		logPrior -= 2.0 * log( alpha[k][j] );

		// the second term of equation 46
		logPrior -= Global::modelBeta * pow( Global::modelGamma, ( double )k ) / alpha[k][j];

	}

	return logPrior;
}

double ModelLearning::calc_lposterior_alphas( double** alpha, size_t k ){
	// calculate partial log posterior of alphas for the order k, due to equation 50 in the theory

	double logPosterior = 0.0;
	float*** v = motif_->getV();
	float** v_bg = bg_->getV();
	size_t N = posSeqs_.size() - 1;

	for( size_t j = 0; j < W_; j++ ){

		// the prior
		logPosterior -= 2.0 * log( alpha[k][j] );

		// the prior
		logPosterior -= Global::modelBeta * pow( Global::modelGamma, ( double )k ) / alpha[k][j];

		//
		if( k == 0 ){

			for( size_t y = 0; y < Y_[k]; y++ ){

				// the third term of equation 50
				logPosterior += boost::math::lgamma( alpha[k][j] );

				// the first part of the forth term of equation 50
				logPosterior += boost::math::lgamma( static_cast<double>( n_[k][y][j] ) + alpha[k][j] * v_bg[k][y] );

				// the second part of the forth term of equation 50
				logPosterior -= boost::math::lgamma( alpha[k][j] * v_bg[k][y] );

			}

			// the fifth term of equation 50
			logPosterior -= boost::math::lgamma( static_cast<double>( N ) + alpha[k][j] );

		} else {

			// the third, forth and fifth terms of equation 50
			for( size_t y = 0; y < Y_[k]; y++ ){

				// !!! important: correction for the occasions when zero or one k-mer is present
				if( n_[k-1][y][j-1] <= 1 ){

					logPosterior = - 2.0 * log( alpha[k][j] ) - Global::modelBeta * pow( Global::modelGamma, static_cast<double>( k ) ) / alpha[k][j];

				} else {
					// the third term of equation 50
					logPosterior += boost::math::lgamma( alpha[k][j] );

					// the forth term of equation 50
					for( size_t A = 0; A < Y_[1]; A++ ){

						size_t ya = y * Y_[1] + A;

						size_t y2 = ya % Y_[k];

						if( n_[k][ya][j] > 0 ){

							// the first part of the forth term
							logPosterior += boost::math::lgamma( static_cast<double>( n_[k][ya][j] ) + alpha[k][j] * v[k-1][y2][j] );

							// the second part of the forth term
							logPosterior -= boost::math::lgamma( alpha[k][j] * v[k-1][y2][j] );
						}

					}

					// the fifth term
					logPosterior -= boost::math::lgamma( static_cast<double>( n_[k-1][y][j-1] ) + alpha[k][j] );
				}
			}
		}
	}
	return logPosterior;
}

double ModelLearning::calc_llikelihood_alphas( double** alpha, size_t k ){
	// calculate partial log likelihood of alphas due to equation 50 in the theory

	double logLikelihood = 0.0;
	float*** v = motif_->getV();
	float** v_bg = bg_->getV();
	double N = static_cast<double>( posSeqs_.size() ) - 1.0;

	if( k == 0 ){

		for( size_t j = 0; j < W_; j++ ){

			for( size_t y = 0; y < Y_[k]; y++ ){

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

		for( size_t j = 0; j < W_; j++ ){
			// the third, forth and fifth terms of equation 50
			for( size_t y = 0; y < Y_[k]; y++ ){

				// the third term of equation 50
				logLikelihood += boost::math::lgamma( alpha[k][j] );

				// the forth term of equation 50
				for( size_t A = 0; A < Y_[1]; A++ ){

					size_t ya = y * Y_[1] + A;

					size_t y2 = ya % Y_[k];

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

void ModelLearning::write( size_t N ){

	/**
	 * 	 * save BaMM (hyper-)parameters in four flat files:
	 * (1) posSequenceBasename.counts:			refined fractional counts of (k+1)-mers
	 * (2) posSequenceBasename.weights: 		responsibilities, posterior distributions
	 * (3) posSequenceBasename.alphas:			optimized hyper-parameter alphas
	 * (4) posSequenceBasename.positions:		position of motif(s) on each sequence
	 * additional for checking:
	 * (5) posSequenceBasename.lpos:			log posterior of alphas with different orders
	 */

	std::string opath = std::string( Global::outputDirectory ) + '/'
						+ Global::posSequenceBasename + "_motif_" + std::to_string( N );

	// output (k+1)-mer counts n[k][y][j]
	std::string opath_n = opath + ".counts";
	std::ofstream ofile_n( opath_n.c_str() );
	if( Global::EM ){
		for( size_t j = 0; j < W_; j++ ){
			for( size_t k = 0; k < K_+1; k++ ){
				for( size_t y = 0; y < Y_[k+1]; y++ ){
					ofile_n << ( int )n_[k][y][j] << '\t';
				}
				ofile_n << std::endl;
			}
			ofile_n << std::endl;
		}
	} else if( Global::CGS ){
		for( size_t j = 0; j < W_; j++ ){
			for( size_t k = 0; k < K_+1; k++ ){
				for( size_t y = 0; y < Y_[k+1]; y++ ){
					ofile_n << n_z_[k][y][j] << '\t';
				}
				ofile_n << std::endl;
			}
			ofile_n << std::endl;
		}
	}

	// output position(s) of motif(s): pos_[n][i]
	std::string opath_pos = opath + ".positions";
	std::ofstream ofile_pos( opath_pos.c_str() );

	ofile_pos << "seq" << '\t' << "position" << std::endl;

	if( Global::EM ){
		float cutoff = 0.3f;	// threshold for having a motif on the sequence in terms of responsibilities
		for( size_t n = 0; n < posSeqs_.size(); n++ ){

			ofile_pos << n+1 << '\t';

			size_t LW1 = posSeqs_[n]->getL() - W_ + 1;

			for( size_t i = LW1; i > 0; i-- ){
				if( r_[n][i] >= cutoff ){
					ofile_pos << LW1-i << '\t';
				}
			}

			ofile_pos << std::endl;
		}
	} else if( Global::CGS ){
		for( size_t n = 0; n < posSeqs_.size(); n++ ){
			ofile_pos << n+1 << '\t' << z_[n]+1 <<'\t' << z_[n]+W_ << std::endl;
		}
	}


	// output hyper-parameter alphas alpha[k][j]
	std::string opath_alpha = opath + ".alphas";
	std::ofstream ofile_alpha( opath_alpha.c_str() );
	for( size_t k = 0; k < K_+1; k++ ){
		ofile_alpha << "> k=" << k << std::endl;
		for( size_t j = 0; j < W_; j++ ){
			ofile_alpha << std::setprecision( 3 ) << alpha_[k][j] << '\t';
		}
		ofile_alpha << std::endl;
	}

/*	// output responsibilities r[n][i]
	std::string opath_r = opath + ".weights";
	std::ofstream ofile_r( opath_r.c_str() );
	for( size_t n = 0; n < posSeqs_.size(); n++ ){
		ofile_r << std::scientific << std::setprecision( 2 ) << r_[n][0] << ' ';
		int LW1 = posSeqs_[n]->getL() - W_ + 1;
		for( i = LW1; i > 0; i-- ){
			ofile_r << std::setprecision( 2 ) << r_[n][i] << ' ';
		}
		ofile_r << std::endl;
	}*/

	// output log posterior alphas
	std::string opath_lpos = opath + ".lpos";
	std::ofstream ofile_lpos( opath_lpos.c_str() );
	for( size_t k = 0; k < K_+1; k++ ){
		ofile_lpos << std::scientific << std::setprecision( 6 ) << calc_lposterior_alphas( alpha_, k ) << '\t';
	}
}
