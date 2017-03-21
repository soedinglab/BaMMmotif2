/*
 * ModelLearning.cpp
 *
 *  Created on: Dec 2, 2016
 *      Author: wanwan
 */

#include "ModelLearning.h"
#include "SeqGenerator.h"

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
	s_ = ( float** )calloc( Y_[K+1]+1, sizeof( float* ) ) + 1;
	for( y = -1; y < Y_[K+1]; y++ ){
		s_[y] = ( float* )calloc( W, sizeof( float ) );
	}

	// allocate memory for r_[n][i], pos_[n][i], z_[n], normFactor_[n] and initialize them
	r_ = ( float** )calloc( posSeqs_.size(), sizeof( float* ) );
	pos_ = ( float** )calloc( posSeqs_.size(), sizeof( float* ) );
	z_ = ( int* )calloc( posSeqs_.size(), sizeof( int ) );
	normFactor_ = ( float* )calloc( posSeqs_.size(), sizeof( float ) );
	for( size_t n = 0; n < posSeqs_.size(); n++ ){
		LW2 = posSeqs_[n]->getL() - W + 2;
		z_[n] = rand() % LW2;
		r_[n] = ( float* )calloc( LW2, sizeof( float ) );
		pos_[n] = ( float* )calloc( LW2, sizeof( float ) );
	}

	// allocate memory for n_[k][y][j] and probs_[k][y][j] and initialize them
	n_ = ( float*** )calloc( K+1, sizeof( float** ) );
	n_z_ = ( int*** )calloc( K+1, sizeof( int** ) );
	for( k = 0; k < K+1; k++ ){
		n_[k] = ( float** )calloc( Y_[k+1]+1, sizeof( float* ) )+1;
		n_z_[k] = ( int** )calloc( Y_[k+1]+1, sizeof( int* ) )+1;
		for( y = -1; y < Y_[k+1]; y++ ){
			n_[k][y] = ( float* )calloc( W, sizeof( float ) );
			n_z_[k][y] = ( int* )calloc( W, sizeof( int ) );
		}
	}

	// allocate memory for alpha_[k][j] and initialize it
	alpha_ = ( float** )malloc( ( K+1 ) * sizeof( float* ) );
	m1_t_ = ( float** )calloc( K+1, sizeof( float* ) );
	m2_t_ = ( float** )calloc( K+1, sizeof( float* ) );
	prob_a_ = ( float** )calloc( K+1, sizeof( float* ) );
	for( k = 0; k < K+1; k++ ){
		alpha_[k] = ( float* )malloc( W * sizeof( float ) );
		m1_t_[k] = ( float* )calloc( W, sizeof( float ) );
		m2_t_[k] = ( float* )calloc( W, sizeof( float ) );
		prob_a_[k] = ( float* )calloc( W, sizeof( float ) );
		for( j = 0; j < W; j++ ){
			alpha_[k][j] = Global::modelAlpha[k];
		}
	}

}

ModelLearning::~ModelLearning(){

	int K = Global::modelOrder;
	for( int y = -1; y < Y_[K+1]; y++ ){
		free( s_[y] );
	}
	free( s_-1 );

	for( size_t n = 0; n < posSeqs_.size(); n++ ){
		free( r_[n] );
		free( pos_[n] );
	}
	free( r_ );
	free( pos_ );
	free( z_ );
	free( normFactor_ );

	for( int k = 0; k < K+1; k++ ){
		for( int y = -1; y < Y_[k+1]; y++ ){
			free( n_[k][y] );
			free( n_z_[k][y] );
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
	llikelihood_ = 0.0f;										// reset log likelihood

	for( size_t n = 0; n < posSeqs_.size(); n++ ){
		normFactor_[n] = 0.0f;
	}
	// calculate responsibilities r_[n][i] at position i in sequence n
	// n runs over all sequences
	// and parallel the code
	#pragma omp parallel for

	for( size_t n = 0; n < posSeqs_.size(); n++ ){

		int L = posSeqs_[n]->getL();
		int LW1 = L - W + 1;
		int LW2 = L - W + 2;
		int* kmer = posSeqs_[n]->getKmer();

		// reset r_[n][i] and pos_[n][i]
		for( int i = 0; i < LW2; i++ ){
			r_[n][i] = 1.0f;
			pos_[n][i] = q_ / static_cast<float>( LW1 );		// p(z_n = i), i > 0
		}
		pos_[n][0] = 1 - q_;

		// when p(z_n > 0)
		// ij = i+j runs over all positions in sequence
		for( int ij = 0; ij < L; ij++ ){

			// extract (K+1)-mer y from positions (i-k,...,i)
			int y = ( kmer[ij] >= 0 ) ? kmer[ij] % Y_[K+1] : -1;

			// j runs over all motif positions
			for( int j = ( 0 > ( ij-L+W ) ? 0 : ij-L+W ); j < ( W < (ij+1) ? W : ij+1 ); j++ ){
				r_[n][L-W+1-ij+j] *= s_[y][j];
			}
		}

		// calculate complete responsibilities and sum them up
		for( int i = 0; i < LW2; i++ ){
			if( r_[n][i] != 0.0f ){
				r_[n][i] *= pos_[n][i];
			}
			normFactor_[n] += r_[n][i];
		}

		// normalize responsibilities
		for( int i = 0; i < LW2; i++ ){
			r_[n][i] /= normFactor_[n];
		}
	}

	for( size_t n = 0; n < posSeqs_.size(); n++ ){
		llikelihood_ += logf( normFactor_[n] );
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
		int* kmer = posSeqs_[n]->getKmer();

		// ij = i+j runs over all positions in x
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
	int K_bg = ( Global::bgModelOrder < K ) ? Global::bgModelOrder : K;
	int W = motif_->getW();
	float eta = Global::eta;						// learning rate for alpha learning
	int k, y, j, n;
	float** alpha_avg = ( float** )calloc( K+1, sizeof( float* ) );
	for( k = 0; k < K+1; k++ ){
		alpha_avg[k] = ( float* )calloc( W, sizeof( float ) );
	}

	// ToDo: write down log posterior and model v after each alpha updating step
/*	std::string opath = std::string( Global::outputDirectory ) + '/'
						+ Global::posSequenceBasename;
	std::string opath_log = opath + ".logposterior";
	std::string opath_condProb = opath + ".condProb";
	std::string opath_alphas = opath + ".alphas";
	std::string opath_gradient = opath + ".gradient";
	std::string opath_vdiff = opath + ".vdiff";
	std::string opath_logAs = opath + ".logAs";

	std::ofstream ofile_log( opath_log.c_str() );
	std::ofstream ofile_vdiff( opath_vdiff.c_str() );
	std::ofstream ofile_condProb( opath_condProb.c_str() );
	std::ofstream ofile_alphas( opath_alphas.c_str() );
	std::ofstream ofile_gradient( opath_gradient.c_str() );
	std::ofstream ofile_logAs( opath_logAs.c_str() );*/

	float*** v_prev = ( float*** )calloc( K+1, sizeof( float** ) );
	for( k = 0; k < K+1; k++ ){
		v_prev[k] = ( float** )calloc( Y_[k+1], sizeof( float* ) );
		for( y = 0; y < Y_[k+1]; y++ ){
			v_prev[k][y] = ( float* )calloc( W, sizeof( float ) );
		}
	}

/*
	// initialize the model with one EStep
	// compute log odd scores s[y][j], log likelihoods of the highest order K
	for( y = 0; y < Y_[K+1]; y++ ){
		for( j = 0; j < W; j++ ){
			int y_bg = y % Y_[K_bg+1];
			s_[y][j] = motif_->getV()[K][y][j] / bg_->getV()[K_bg][y_bg];
		}
	}
	// E-step: calculate posterior
	EStep();
	// extract initial z from the indices of the largest responsibilities
	for( n = 0; n < ( int )posSeqs_.size(); n++ ){
		int LW2 = posSeqs_[n]->getL() - motif_->getW() + 2;
		for( int i = 1; i < LW2; i++ ){
			z_[n] = ( r_[n][i] > r_[n][i-1] ) ? i : i-1;
		}
	}
*/

	// count the k-mers
	// 1. reset n_z_[k][y][j]
	for( k = 0; k < K+1; k++ ){
		for( y = 0; y < Y_[k+1]; y++ ){
			for( j = 0; j < W; j++ ){
				n_z_[k][y][j] = 0;
			}
		}
	}
	// 2. count k-mers for all the orders
	// note: the last sequence is not counted, since during the first iteration
	// of iteration, it is added back.
	for( n = 0; n < ( int )posSeqs_.size()-1; n++ ){
		int* kmer = posSeqs_[n]->getKmer();
		if( z_[n] > 0 ){
			for( k = 0; k < K+1; k++ ){
				for( j = 0; j < W; j++ ){
					y = ( kmer[z_[n]-1+j] >= 0 ) ? kmer[z_[n]-1+j] % Y_[k+1] : -1;
					n_z_[k][y][j]++;
				}
			}
		}
	}

	// iterate over
	while( iterate && iteration < Global::maxCGSIterations ){

		iteration++;

		// get the previous v before updating alphas and z and q
		for( k = 0; k < K+1; k++ ){
			for( y = 0; y < Y_[k+1]; y++ ){
				for( j = 0; j < W; j++ ){
					v_prev[k][y][j] = motif_->getV()[k][y][j];
				}
			}
		}

		// sampling z and q
		if( !Global::noZQSampling )		Gibbs_sampling_z_q();

		// update alphas:
		if( !Global::noAlphaOptimization ){

			// update alphas by stochastic optimization
			stochastic_optimize_alphas( K, W, eta, iteration );

			// update model parameter v
			motif_->updateVz_n( n_z_, alpha_, K );

/*
			// todo: write out the difference of v's for comparison
			for( k = 0; k < K+1; k++ ){
				float v_diff = 0.0f;
				for( y = 0; y < Y_[k+1]; y++ ){
					for( j = 0; j < W; j++ ){
						v_diff += fabsf( motif_->getV()[k][y][j] - v_prev[k][y][j] );
					}
				}
				ofile_vdiff << v_diff << '\t';
			}
			ofile_vdiff << std::endl;

			for( j = 0; j < W; j++ ){
				ofile_alphas << alpha_[K][j] << '\t';
				ofile_gradient << calc_gradient_alphas( alpha_, K, j ) << '\t';
				ofile_logAs << calc_logCondProb_a( alpha_, K, j ) << '\t';

			}
			ofile_alphas << std::endl;
			ofile_gradient << std::endl;
			ofile_logAs << std::endl;

			ofile_condProb << expf( calc_logCondProb_alphas( alpha_, K ) ) << '\t';
			ofile_log << calc_logCondProb_alphas( alpha_, K ) << '\t'
					<< calc_prior_alphas( alpha_, K ) << '\t'
					<< calc_likelihood_alphas( alpha_, K ) << '\t';
			ofile_log << std::endl;
			ofile_condProb << std::endl;
*/

		} else if( Global::alphaSampling ){
			// sample alpha in the exponential space using Metreopolis-Hastings algorithm
			GibbsMH_sampling_alphas();

		} else if( Global::debugAlphas ){

			// debug the optimization of alphas
			if( iteration <= 150 ){
				// update alphas by stochastic optimization
				stochastic_optimize_alphas( K, W, eta, iteration );

				// update model parameter v
				motif_->updateVz_n( n_z_, alpha_, K );

			} else {

				// set up a sequential numbers for alphas with the highest order
				for( j = 0; j < W; j++ ){
					alpha_[K][j] = ( float )iteration - 150.0f;
				}

				// update model parameter v
				motif_->updateVz_n( n_z_, alpha_, K );

/*				// todo: write out the difference of v's for comparison
				for( k = 0; k < K+1; k++ ){
					float v_diff = 0.0f;
					for( y = 0; y < Y_[k+1]; y++ ){
						for( j = 0; j < W; j++ ){
							v_diff += fabsf( motif_->getV()[k][y][j] - v_prev[k][y][j] );
						}
					}
					ofile_vdiff << v_diff << '\t';
				}
				ofile_vdiff << std::endl;

				for( j = 0; j < W; j++ ){
					ofile_alphas << alpha_[K][j] << '\t';
					ofile_gradient << calc_gradient_alphas( alpha_, K, j ) << '\t';
					ofile_logAs << calc_logCondProb_a( alpha_, K, j ) << '\t';

				}
				ofile_alphas << std::endl;
				ofile_gradient << std::endl;
				ofile_logAs << std::endl;

				ofile_condProb << expf( calc_logCondProb_alphas( alpha_, K ) ) << '\t';
				ofile_log << calc_logCondProb_alphas( alpha_, K ) << '\t'
						<< calc_prior_alphas( alpha_, K ) << '\t'
						<< calc_likelihood_alphas( alpha_, K ) << '\t';
				ofile_log << std::endl;
				ofile_condProb << std::endl;*/
			}

		}

		// only for writing out model after each iteration:
/*		motif_->calculateP();
		motif_->write( CGSIterations_ );*/

		// get the sum of alpha[k][j] at the last five iterations
		if( iteration > Global::maxCGSIterations - 5 ){
			for( k = 0; k < K+1; k++ ){
				for( j = 0; j < W; j++ ){
					alpha_avg[k][j] += alpha_[k][j];
				}
			}
		}

	}

/*	// obtaining a motif model
	// get the average alpha[k][j] from the last five iterations
	for( k = 0; k < K+1; k++ ){
		for( j = 0; j < W; j++ ){
			alpha_avg[k][j] /= 5.0f;
		}
	}
	// run five steps of EM to optimize the final model with
	// the optimum model parameters v's and the fixed alphas
	for( size_t step = 0; step < 3; step++ ){
		// compute log odd scores s[y][j], log likelihoods of the highest order K
		for( y = 0; y < Y_[K+1]; y++ ){
			for( j = 0; j < W; j++ ){
				int y_bg = y % Y_[K_bg+1];
				s_[y][j] = motif_->getV()[K][y][j] / bg_->getV()[K_bg][y_bg];
			}
		}

		// E-step: calculate posterior
		EStep();

		// M-step: update model parameters
		MStep();
	}*/


	// calculate probabilities
	motif_->calculateP();

	fprintf( stdout, "\n--- Runtime for Collapsed Gibbs sampling: %.4f seconds ---\n", ( ( float )( clock() - t0 ) ) / CLOCKS_PER_SEC );

	// free the memory:
	for( k = 0; k < K+1; k++ ){
		for( y = 0; y < Y_[k+1]; y++ ){
			free( v_prev[k][y] );
		}
		free( v_prev[k] );
		free( alpha_avg[k] );
	}
	free( v_prev );
	free( alpha_avg );

}

void ModelLearning::Gibbs_sampling_z_q(){

	int N = ( int )posSeqs_.size();
	int W = motif_->getW();
	int K = Global::modelOrder;
	int K_bg = ( Global::bgModelOrder < K ) ? Global::bgModelOrder : K;
	int k, y, y_prev, y_bg, j, i, LW1, n, n_prev;
	int N_0 = 0;								// counts of sequences which do not contain motifs.

	float** v_bg = bg_->getV();

	// sampling z:
	// loop over all sequences and drop one sequence each time and update r
	for( n = 0; n < N; n++ ){

		LW1 = posSeqs_[n]->getL() - W + 1;
		int* kmer = posSeqs_[n]->getKmer();

		// calculate positional prior:
		pos_[n][0] = 1.0f - q_;
		for( i = 1; i <= LW1; i++ ){
			pos_[n][i] = q_ / ( float )LW1;
		}

		// count k-mers at position z_i+j except the n'th sequence
		// remove the k-mer counts from the current sequence with old z
		if( z_[n] > 0 ){
			for( k = 0; k < K+1; k++ ){
				for( j = 0; j < W; j++ ){
					y = ( kmer[z_[n]-1+j] >= 0 ) ? kmer[z_[n]-1+j] % Y_[k+1] : -1;
					n_z_[k][y][j]--;
				}
			}
		}

		// get the index of the previous sequence
		n_prev = ( n > 0 ) ? n-1 : N-1;

		// add the k-mer counts from the previous sequence with updated z
		int* kmer_prev = posSeqs_[n_prev]->getKmer();
		if( z_[n_prev] > 0 ){
			for( k = 0; k < K+1; k++ ){
				for( j = 0; j < W; j++ ){
					y_prev = ( kmer_prev[z_[n_prev]-1+j] >= 0 ) ? kmer_prev[z_[n_prev]-1+j] % Y_[k+1] : -1;
					n_z_[k][y_prev][j]++;
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
		// todo: could be parallelized by extracting 8 sequences at once
		for( i = 1; i <= LW1; i++ ){
			r_[n][i] = 1.0f;
			for( j = 0; j < W; j++ ){
				// extract k-mers on the motif at position i over W of the n'th sequence
				y = ( kmer[i-1+j] >= 0 ) ? kmer[i-1+j] % Y_[K+1] : -1;
				r_[n][i] *= s_[y][j];
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

		// draw a sample z randomly
		z_[n] = posterior_dist( Global::rngx );

		if( z_[n] == 0 ) N_0++;
	}

	// sampling q:
	// draw two random numbers Q and P from Gamma distribution
	std::gamma_distribution<> P_Gamma_dist( N_0 + 1, 1 );
	std::gamma_distribution<> Q_Gamma_dist( N - N_0 + 1, 1 );
	// draw a sample for P
	double P = P_Gamma_dist( Global::rngx );
	// draw a sample for Q
	double Q = Q_Gamma_dist( Global::rngx );
	// calculate q_
	q_ = ( float )Q / ( float )( Q + P );

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
			// Note: here change the sign in front of eta from '-' to '+'
			a += eta * m1 / ( ( sqrtf( m2 ) + epsilon ) * sqrtf( t ) );

			alpha_[k][j] = expf( a );
		}

	}

}

void ModelLearning::GibbsMH_sampling_alphas(){
	// Gibbs sampling alphas in exponential space with Metropolis-Hastings algorithm

	int K = Global::modelOrder;
	int W = motif_->getW();

	for( int k = 0; k < K+1; k++ ){

		for( int j = 0; j < W; j++ ){

			// draw 10 times in a row and take record of the last accepted sample
			for( int step = 0; step < 10; step++ ){

				float lprob_a_prev = calc_logCondProb_a( alpha_, k, j );

				// draw a new a from N(a, 1)
				std::normal_distribution<float> norm_dist( logf( alpha_[k][j] ), 1.0f );
				alpha_[k][j] = expf( norm_dist( Global::rngx ) );

				// accept this trial sample with a probability
				float lprob_a_new = calc_logCondProb_a( alpha_, k, j );
				prob_a_[k][j] = ( lprob_a_new < lprob_a_prev ) ? expf( lprob_a_new - lprob_a_prev ) : 1;
			}
		}
	}

}

float ModelLearning::calc_gradient_alphas( float** alpha, int k, int j ){
	// calculate partial gradient of the log posterior of alphas due to equation 47 in the theory
	// Note that j >= k

	float gradient = 0.0f;
	float*** v = motif_->getV();
	float** v_bg = bg_->getV();
	int N = static_cast<int>( posSeqs_.size() ) - 1;
	int y, y2;


	// the first term of equation 47
	gradient -= 2.0f / alpha[k][j];

	// the second term of equation 47
	gradient += Global::modelBeta * powf( Global::modelGamma, ( float )k ) / powf( alpha[k][j], 2.f );

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

		if( k == 0 ){

			gradient -= boost::math::digamma( ( float )N + alpha[k][j] );

		} else if( j == 0 ){

			gradient -= boost::math::digamma( ( float )n_z_[k-1][y][j] / ( float )Y_[1] + alpha[k][j] );

		} else {

			gradient -= boost::math::digamma( ( float )n_z_[k-1][y][j-1] + alpha[k][j] );

		}
	}
	return gradient;
}

float ModelLearning::calc_logCondProb_alphas( float** alpha, int k ){
	// calculate partial log conditional probabilities of alphas due to equation 46 in the theory

	float logCondProb = 0.0f;
	float*** v = motif_->getV();
	float** v_bg = bg_->getV();
	int W = motif_->getW();
	int N = static_cast<int>( posSeqs_.size() ) - 1;
	int j, y, y2;

	for( j = 0; j < W; j++ ){

		// the first term of equation 46
		logCondProb -= 2.0f * logf( alpha[k][j] );

		// the second term of equation 46
		logCondProb -= Global::modelBeta * powf( Global::modelGamma, ( float )k ) / alpha[k][j];

		// the third term of equation 46
		logCondProb += ( float )ipow( Y_[1], k ) * boost::math::lgamma( alpha[k][j] );

		// the forth term of equation 46
		for( y = 0; y < Y_[k+1]; y++ ){
			y2 = y % Y_[k];									// cutoff the first nucleotide in the (k+1)-mer
			if( k == 0 ){
				// the first part
				logCondProb += boost::math::lgamma( ( float )n_z_[k][y][j] + alpha[k][j] * v_bg[k][y] );

				// the second part
				logCondProb -= boost::math::lgamma( alpha[k][j] * v_bg[k][y] );


			} else {
				// the first part
				logCondProb += boost::math::lgamma( ( float )n_z_[k][y][j] + alpha[k][j] * v[k-1][y2][j] );

				// the second part
				logCondProb -= boost::math::lgamma( alpha[k][j] * v[k-1][y2][j] );

			}
		}

		// the last term of equation 46
		for( y = 0; y < Y_[k]; y++ ){

			if( j == 0 || k == 0 ){

				logCondProb -= boost::math::lgamma( ( float )N / ( float )Y_[k] + alpha[k][j] );

			} else {

				logCondProb -= boost::math::lgamma( ( float )n_z_[k-1][y][j-1] + alpha[k][j] );

			}
		}
	}

	return logCondProb;
}

float ModelLearning::calc_logCondProb_a( float** alpha, int k, int j ){
	// calculate partial log conditional probabilities of a's due to equation 50 in the theory

	float logCondProbA = 0.0f;
	float*** v = motif_->getV();
	float** v_bg = bg_->getV();
	int N = static_cast<int>( posSeqs_.size() ) - 1;
	int y, y2;


	// the first term of equation 50
	logCondProbA -= logf( alpha[k][j] );

	// the second term of equation 50
	logCondProbA -= Global::modelBeta * powf( Global::modelGamma, ( float )k ) / alpha[k][j];

	// the third term of equation 50
	logCondProbA += ( float )ipow( Y_[1], k ) * boost::math::lgamma( alpha[k][j] );

	// the forth term of equation 50
	for( y = 0; y < Y_[k+1]; y++ ){

		// cut off the first nucleotide in the (k+1)-mer
		y2 = y % Y_[k];

		if( k == 0 ){
			// the first part
			logCondProbA += boost::math::lgamma( ( float )n_z_[k][y][j] + alpha[k][j] * v_bg[k][y] );

			// the second part
			logCondProbA -= boost::math::lgamma( alpha[k][j] * v_bg[k][y] );

		} else {
			// the first part
			logCondProbA += boost::math::lgamma( ( float )n_z_[k][y][j] + alpha[k][j] * v[k-1][y2][j] );

			// the second part
			logCondProbA -= boost::math::lgamma( alpha[k][j] * v[k-1][y2][j] );

		}
	}

	// the last term of equation 46
	for( y = 0; y < Y_[k]; y++ ){

		if( j == 0 || k == 0 ){

			logCondProbA -= boost::math::lgamma( ( float )N + alpha[k][j] );

		} else {

			logCondProbA -= boost::math::lgamma( ( float )n_z_[k-1][y][j-1] + alpha[k][j] );

		}
	}

	return logCondProbA;
}

float ModelLearning::calc_prior_alphas( float** alpha, int k ){
	// calculate partial log conditional probabilities of alphas due to equation 46 in the theory

	float logPrior = 0.0f;
	int W = motif_->getW();

	for( int j = 0; j < W; j++ ){

		// the first term of equation 46
		logPrior -= 2.0f * logf( alpha[k][j] );

		// the second term of equation 46
		logPrior -= Global::modelBeta * powf( Global::modelGamma, ( float )k ) / alpha[k][j];

	}

	return logPrior;
}

float ModelLearning::calc_likelihood_alphas( float** alpha, int k ){
	// calculate partial log likelihood of alphas due to equation 46 in the theory

	float logLikelihood = 0.0f;
	float*** v = motif_->getV();
	float** v_bg = bg_->getV();
	float N = static_cast<float>( posSeqs_.size() ) - 1.0f;
	int W = motif_->getW();
	int y, y2, j;

	for( j = 0; j < W; j++ ){

		// the third term of equation 46
		logLikelihood += ( float )ipow( Y_[1], k ) * boost::math::lgamma( alpha[k][j] );

		// the forth term of equation 46
		for( y = 0; y < Y_[k+1]; y++ ){

			y2 = y % Y_[k];									// cutoff the first nucleotide in the (k+1)-mer

			if( k == 0 ){
				// the first part
				logLikelihood += boost::math::lgamma( ( float )n_z_[k][y][j] + alpha[k][j] * v_bg[k][y] );

				// the second part
				logLikelihood -= boost::math::lgamma( alpha[k][j] * v_bg[k][y] );

			} else {
				// the first part
				logLikelihood += boost::math::lgamma( ( float )n_z_[k][y][j] + alpha[k][j] * v[k-1][y2][j] );

				// the second part
				logLikelihood -= boost::math::lgamma( alpha[k][j] * v[k-1][y2][j] );

			}

			// the missing part for v_bg
//			logLikelihood -= ( float )n_z_[k][y][j] * logf( v_bg[k][y] );

		}

		// the last term of equation 46
		for( y = 0; y < Y_[k]; y++ ){
			if( j == 0 || k == 0 ){
				logLikelihood -= boost::math::lgamma( N / ( float )Y_[k] + alpha[k][j] );
			} else {
				logLikelihood -= boost::math::lgamma( ( float )n_z_[k-1][y][j-1] + alpha[k][j] );
			}
		}

	}

	return logLikelihood;
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
			std::cout << "alpha = " << alphas[k][j] << '\t';
			// calculate the delta of log posterior
			alphas[k][j] += stepsize;
			diff = calc_logCondProb_alphas( alphas, k );

			alphas[k][j] -= 2.0f * stepsize;
			diff -= calc_logCondProb_alphas( alphas, k );

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

float ModelLearning::debug_first_term_of_derivative_with_prior( int k, int j ){
	// interpretation of the partial derivative due to e.q. 49
	// sum up prior and the first and third term of derivative
	// Note that k > 1

	float*** v = motif_->getV();

	float N = ( float )posSeqs_.size();

	float s = 0.0f;

	float third_term = 0.0f;

	// calculate the prior
	float prior = -2.0f / alpha_[k][j] + Global::modelBeta *
			powf( Global::modelGamma, ( float )k ) / powf( alpha_[k][j], 2.0f );

	for( int y = 0; y < Y_[k+1]; y++ ){

		int y2 = y % Y_[k];

		int yk = y / Y_[1];

		s -= v[k-1][j][y2] * ( logf( v[k][y][j] / v[k-1][y2][j] ) );

		if( j > 0 ){
			third_term = 0.5f * ( v[k][y][j] - v[k-1][y2][j] ) / ( v[k][y][j] * ( ( float )n_z_[k-1][yk][j-1] + 1.0e-8f ) );

		} else {
			third_term = 0.5f * ( v[k][y][j] - v[k-1][y2][j] ) / ( v[k][y][j] * N );
		}
	}

	s -= prior;

	s -= third_term;

	return s;

}

float ModelLearning::debug_second_term_of_derivative( int k, int j ){
	// interpretation of the partial derivative due to e.q. 49
	// calculate the second term of derivative
	// Note that k > 1

	float N = ( float )posSeqs_.size();

	float s = 0.0f;

	for( int y = 0; y < Y_[k]; y++ ){
		if( j > 0 ){
			s += 1.5f * ( float )n_z_[k-1][y][j-1] / ( ( ( float )n_z_[k-1][y][j-1] + alpha_[k][j] ) * alpha_[k][j] );
		} else {
			s += 1.5f * N / ( ( N + alpha_[k][j] ) * alpha_[k][j] );
		}
	}

	return s;
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
//		std::string opath_r = opath + ".EMweight";
		std::string opath_pos = opath + ".EMposition";
//		std::ofstream ofile_r( opath_r.c_str() );
		std::ofstream ofile_pos( opath_pos.c_str() );
		ofile_pos << "seq" << '\t' << "positions" << std::endl;
		float cutoff = 0.3f;									// threshold for having a motif on the sequence in term of responsibilities
		for( size_t n = 0; n < posSeqs_.size(); n++ ){
			int LW1 = posSeqs_[n]->getL() - W + 1;
//			ofile_r << std::scientific << std::setprecision( 2 ) << r_[n][0] << ' ';	// print out the responsibility of not having a motif on the sequence
			ofile_pos << n+1 << '\t';							// print out the sequence number
			for( i = LW1; i > 0; i-- ){
//				ofile_r << std::setprecision( 2 ) << r_[n][i] << ' ';
				if( r_[n][i] >= cutoff ){
					ofile_pos << LW1-i << '\t';
				}
			}
//			ofile_r << std::endl;
			ofile_pos << std::endl;
		}

/*
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
*/

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

/*		// output responsibilities r[n][i]
		std::string opath_r = opath + ".CGSweight";
		std::ofstream ofile_r( opath_r.c_str() );
		for( size_t n = 0; n < posSeqs_.size(); n++ ){
			for( i = 0; i < posSeqs_[n]->getL()-W+2; i++ ){
				ofile_r << std::scientific << std::setprecision( 2 ) << r_[n][i] << ' ';
			}
			ofile_r << std::endl;
		}*/

		// output parameter alphas alpha[k][j]
		std::string opath_alpha = opath + ".CGSalpha";
		std::ofstream ofile_alpha( opath_alpha.c_str() );
		for( k = 0; k < K+1; k++ ){
			ofile_alpha << "k = " << k << std::endl;
			for( j = 0; j < W; j++ ){
				ofile_alpha << std::setprecision( 3 ) << alpha_[k][j] << '\t';
			}
			ofile_alpha << std::endl;
		}

		// output positions of motifs z_[n]
		std::string opath_z = opath + ".CGSposition";
		std::ofstream ofile_z( opath_z.c_str() );
		ofile_z << "seq" << '\t' << "start" <<'\t' << "end" << std::endl;
		for( size_t n = 0; n < posSeqs_.size(); n++ ){
			ofile_z << n+1 << '\t' << z_[n]+1 <<'\t' << z_[n]+W << std::endl;
		}
	}
}
