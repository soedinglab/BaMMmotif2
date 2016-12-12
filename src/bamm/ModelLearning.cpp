/*
 * ModelLearning.cpp
 *
 *  Created on: Dec 2, 2016
 *      Author: wanwan
 */

#include "ModelLearning.h"

#include <random>
#include <cmath>

ModelLearning::ModelLearning( Motif* motif, BackgroundModel* bg, std::vector<int> folds ){

	motif_ = motif;
	bg_ = bg;

	// define parameters
	int K = Global::modelOrder;
	int K_bg = Global::bgModelOrder;
	int W = motif_->getW();
	int y, k, j, LW2;
	for( k = 0; k < std::max( K+2, K_bg+2 ); k++ ){
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
		z_[n] = rand() % LW2;
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
	alpha_ = ( float** )malloc( ( K+1 ) * sizeof( float* ) );
	for( k = 0; k < K+1; k++ ){
		alpha_[k] = ( float* )malloc( W * sizeof( float ) );
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
	}
	free( n_ );
	free( n_z_ );
	free( alpha_ );

}

int ModelLearning::EMlearning(){

	printf( " ______\n"
			"|      |\n"
			"|  EM  |\n"
			"|______|\n\n" );

	clock_t t0 = clock();
	bool iterate = true;									// flag for iterating before convergence
	int W = motif_->getW();
	int K = Global::modelOrder;
	int K_bg = Global::bgModelOrder;

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
				s_[y][j] = motif_->getV()[K][y][j] / bg_->getV()[std::min( K, K_bg )][y_bg];
			}
		}

		// E-step: calculate posterior
		EM_EStep();

		// M-step: update model parameters
		EM_MStep();

		// * optional: optimize parameter alpha
		if( !Global::noAlphaOptimization )	EM_optimizeAlphas( K, W );

		// * optional: optimize parameter q
		if( !Global::noQOptimization )		EM_optimize_q();

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

void ModelLearning::EM_EStep(){

	int K = Global::modelOrder;
	int W = motif_->getW();
	llikelihood_ = 0.0f;										// reset log likelihood

	// calculate responsibilities r_[n][i] at position i in sequence n
	for( size_t n = 0; n < posSeqs_.size(); n++ ){				// n runs over all sequences
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
		for( int ij = 0; ij < L; ij++ ){						// ij = i+j runs over all positions in sequence
			int y = posSeqs_[n]->extractKmer( ij, std::min( ij, K ) );			// extract (k+1)-mer y from positions (i-k,...,i)
			for( int j = std::max( 0, ij-L+W ); j < std::min( W, ij+1 ); j++ ){	// j runs over all motif positions
				if( y != -1 ){									// skip 'N' and other unknown alphabets
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

void ModelLearning::EM_MStep(){

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
	for( size_t n = 0; n < posSeqs_.size(); n++ ){				// n runs over all sequences
		int L = posSeqs_[n]->getL();
		int LW1 = L - W + 1;
		for( int ij = 0; ij < L; ij++ ){						// ij = i+j runs over all positions in x
			int y = posSeqs_[n]->extractKmer( ij, std::min( ij, K ) );
			for( int j = std::max( 0, ij-L+W ); j < std::min( W, ij+1 ); j++ ){
				if( y != -1 && ( ij-j ) < LW1 ){				// skip 'N' and other unknown alphabets
					n_[K][y][j] += r_[n][L-W+1-ij+j];
				}
			}
		}
	}

	// compute fractional occurrence counts from higher to lower order
	for( int k = K; k > 0; k-- ){								// k runs over all orders
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

void ModelLearning::EM_optimizeAlphas( int K, int W ){
/*
	float optim_alpha = zbrent( *this, &ModelLearning::EMcalcGrad_Qfunc, min_brent, max_brent, tolerance, K );

	// only update in case a root is bracketed
	if( optim_alpha > 0 ){
		for( int j = 0; j < W; j++ ){
			alpha_[K][j] = optim_alpha;
		}
		motif_->updateV( n_, alpha_, K );

	}*/
}

void ModelLearning::EM_optimize_q(){

	// optimize hyper-parameter q
	// motif.updateV()
}

float ModelLearning::EM_calcLogPriors( int K ){

	int j,y,y2;
	float lPriors = 0.0f;
	int A = Alphabet::getSize();
	int W = motif_->getW();
	float*** v_motif = motif_->getV();
	float** v_bg = bg_->getV();

	// the second and third parts of log Posterior Probability
	for( j = 0; j < W; j++ ){
		// the second part
		lPriors += ( float )Y_[K] * lgammaf( alpha_[K][j] + ( float )A );
		// the second and third terms
		for( y = 0; y < Y_[K+1]; y++ ){
			// the second term
			y2 = y % Y_[K];                         				// cut off the first nucleotide in (k+1)-mer y
			if( K == 0 ){
				lPriors -= lgammaf( alpha_[K][j] * v_bg[K][y] + 1.0f );
				// the third term
				lPriors += alpha_[K][j] * v_bg[K][y] * logf( v_motif[K][y][j] );
			}
			if( K > 0 ){
				lPriors -= lgammaf( alpha_[K][j] * v_motif[K-1][y2][j] + 1.0f );
				// the third term
				lPriors += alpha_[K][j] * v_motif[K-1][y2][j] * logf( v_motif[K][y][j] );
			}
		}
		// the forth part
		lPriors += ( - 2.0f * logf( alpha_[K][j] ) - Global::modelBeta * powf( Global::modelGamma, ( float )K ) /
				alpha_[K][j] + logf( Global::modelBeta * powf( Global::modelGamma, ( float )K ) ) );
	}

	return lPriors;
}

float ModelLearning::EM_calcLogPosterior( int K ){

	return llikelihood_ + EM_calcLogPriors( K );
}

float ModelLearning::EM_calcQfunc( int K ){

	int L;
	float prior_i, prior_0 = 1 - q_;
	int W = motif_->getW();
	float*** v_motif = motif_->getV();
	float Qfunc = 0.0f;

    // the likelihood part of Q function
	for( int y = 0; y < Y_[K+1]; y++ ){
		for( int j = 0; j < W; j++ ){
			Qfunc += n_[K][y][j] * logf( v_motif[K][y][j] );
		}
	}

	// the constant t of Q function
	for( size_t n = 0; n < posSeqs_.size(); n++ ){
		L = posSeqs_[n]->getL();
		prior_i = q_ / static_cast<float>( L - W + 1 );
		Qfunc += ( prior_0 * logf( prior_0 ) + q_ * logf( prior_i ) );
	}

	// the priors of Q function
	Qfunc += EM_calcLogPriors( K );

	return Qfunc;
}

float ModelLearning::EM_calcGrad_Qfunc( float alpha ,int K, int W, int A ){

	int j, y, y2;
	float Qfunc_grad = 0.0f;
    float*** v_motif = motif_->getV();
	float** v_bg = bg_->getV();

    if( K == 0 ){
        for( j = 0; j < W; j++ ){
            // first term of gradient
            Qfunc_grad += digammaf( alpha + (float) A );
            for( y = 0; y < Y_[K+1]; y++ ){
                // second term of gradient
                Qfunc_grad += -( v_bg[K][y] * ( digammaf( alpha * v_bg[K][y] + 1.0f ) - logf( v_motif[K][y][j] )));
            }
            //third term of gradient
            Qfunc_grad += - 2.0f / alpha + Global::modelBeta / powf( alpha , 2.0f ) ;
        }
    }
    if( K > 0 ){
        for( j = 0; j < W; j++ ){
            // first term of gradient
            Qfunc_grad += ( float )Y_[K] * digammaf( alpha + (float) A );
            for( y = 0; y < Y_[K+1]; y++ ){
                y2 = y % Y_[K];
                // second term of gradient
                Qfunc_grad += -( v_motif[K-1][y2][j] * ( digammaf( alpha * v_motif[K-1][y2][j] + 1.0f ) - logf( v_motif[K][y][j] )));
            }
            //third term of gradient
            Qfunc_grad += - 2.0f / alpha + ( Global::modelBeta * powf( Global::modelGamma, (float) K ) ) / powf( alpha , 2.0f ) ;
        }

    }
	return Qfunc_grad;
}

void ModelLearning::GibbsSampling(){

	printf( " ___________________________\n"
			"|                           |\n"
			"|  Collapsed Gibbs sampler  |\n"
			"|___________________________|\n\n" );

	clock_t t0 = clock();
	bool iterate = true;							// flag for iterating before convergence

	int K = Global::modelOrder;
	int W = motif_->getW();

	// parameters for alpha learning
	int timestep = 0;								// timestep = iteration

	// iterate over
	while( iterate && timestep < Global::maxCGSIterations ){

		timestep++;

		if( Global::verbose ){
			std::cout << timestep << " iteration:\n";
		}

		// sampling z and q
		CGS_sampling_z_q();

		// only for writing out model after each iteration:
/*		motif_->calculateP();
		motif_->write( CGSIterations_ );*/

		// update model parameter v
		motif_->updateV( n_z_, alpha_, K );

		// optimize hyper-parameter alpha
		if( !Global::noAlphaUpdating )	CGS_updateAlphas( K, W );

		// test alphas
		// a toy example:
		float** toys = new float*[3];
		for( int i = 0; i < 3; i++ ){
			toys[i] = new float[10];
			for( int b = 0; b < 10; b++ ){
				toys[i][b] = 10.0f;
			}
		}
		testAlphaUpdate( toys, 2, 10 );
		delete[] toys;
//		testAlphaUpdate( alpha_, K, W );

	}

	// obtaining a motif model

	// calculate probabilities
	motif_->calculateP();

	fprintf( stdout, "\n--- Runtime for Collapsed Gibbs sampling: %.4f seconds ---\n", ( ( float )( clock() - t0 ) ) / CLOCKS_PER_SEC );
}

void ModelLearning::CGS_sampling_z_q(){

	int N = ( int )posSeqs_.size();
	int W = motif_->getW();
	int K = Global::modelOrder;
	int K_bg = Global::bgModelOrder;
	int k, y, y_prev, y2, y_bg, j, i, LW1, n, n_prev;
	int N_0 = 0;								// counts of sequences which do not contain motifs.

	// reset n_z_[k][y][j]
	for( k = 0; k < K+1; k++ ){
		for( y = 0; y < Y_[k+1]; y++ ){
			for( j = 0; j < W; j++ ){
				n_z_[k][y][j] = 0;
			}
		}
	}
	// count k-mers for the highest order K
	for( i = 0; i < N; i++ ){
		for( j = 0; j < W; j++ ){
			y = posSeqs_[i]->extractKmer( z_[i]-1+j, std::min( z_[i]-1+j, K ) );
			if( y >= 0 ){
				n_z_[K][y][j]++;
			}
		}
	}
	// calculate nz for lower order k
	for( k = K; k > 0; k-- ){					// k runs over all lower orders
		for( y = 0; y < Y_[k+1]; y++ ){
			y2 = y % Y_[k];						// cut off the first nucleotide in (k+1)-mer
			for( j = 0; j < W; j++ ){
				n_z_[k-1][y2][j] += n_z_[k][y][j];
			}
		}
	}

	// sampling z:
	// loop over all sequences and drop one sequence each time and update r
	for( n = 0; n < N; n++ ){

		LW1 = posSeqs_[n]->getL() - W + 1;

		// calculate positional prior:
		pos_[n][0] = 1.0f - q_;
		for( i = 1; i <= LW1; i++ ){
			pos_[n][i] = q_ / ( float )LW1;
		}

		// count k-mers at position z[i]+j except the n'th sequence
		for( k = 0; k < K+1; k++ ){
			for( j = 0; j < W; j++ ){
				if( z_[n] != 0 ){
					// remove the k-mer counts for the current sequence with old z
					y = posSeqs_[n]->extractKmer( z_[n]-1+j, std::min( z_[n]-1+j, k ) );
					if( y >= 0 ){
						n_z_[k][y][j]--;
					}
				}
				if( n > 0 ){
					n_prev = n-1;
					if( z_[n_prev] != 0 ){
						// add the k-mer counts for the previous sequence with updated z
						y_prev = posSeqs_[n-1]->extractKmer( z_[n_prev]-1+j, std::min( z_[n_prev]-1+j, k ) );
						if( y_prev >= 0 ){
							n_z_[k][y_prev][j]++;
						}
					}
				}
			}
		}

		// updated model parameters v excluding the n'th sequence
		motif_->updateV( n_z_, alpha_, K );

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
				y = posSeqs_[n]->extractKmer( i-1+j, std::min( i-1+j, K ) );
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
			posterior_array.push_back( r_[n][i] );
		}

		// draw a new position z from discrete posterior distribution
		std::discrete_distribution<> posterior_dist( posterior_array.begin(), posterior_array.end() );
		std::random_device rand;						// pick a random number
		z_[n] = posterior_dist( rand );					// draw a sample z

		if( z_[n] == 0 ) N_0++;
	}

	// sampling q:
	// draw two random numbers Q and P from Gamma distribution
	std::gamma_distribution<> P_Gamma_dist( N_0 + 1, 1 );
	std::gamma_distribution<> Q_Gamma_dist( N - N_0 + 1, 1 );
	std::random_device rand1;							// pick a random number
	std::random_device rand2;							// pick another random number
	double P = P_Gamma_dist( rand1 );					// draw a sample for P
	double Q = Q_Gamma_dist( rand2 );					// draw a sample for Q

	q_ = ( float ) Q / ( float )( Q + P );				// calculate q_
//	q_ = ( float )( N - N_0 ) / ( float )N;

	if( Global::verbose ){
		// checking z values from the first 20 sequences
//		for( int m = 0; m < 20; m++ ) std::cout << z_[m] << '\t';
//		std::cout << N_0 << " sequences do not have motif. q = " << q_ << "\n";
	}
}

void ModelLearning::CGS_updateAlphas( int K, int W ){

/*
	// update alphas using stochastic optimization algorithm ADAM (DP Kingma & JL Ba 2015)
	int k, j;
	float eta = 0.0001f;									// learning rate for a
	float beta1 = 0.9f;									// exponential decay rate for the moment estimates
	float beta2 = 0.999f;								// exponential decay rate for the moment estimates
	float epsilon = 0.00000001f;						// cutoff
	float** a;											// a = e^alpha
	float** gradient;									// gradient of log posterior of alpha
	float** m1;											// first moment vector (the mean)
	float** m2;											// second moment vector (the uncentered variance)

	a = ( float** )calloc( K+1, sizeof( float* ) );
	gradient = ( float** )calloc( K+1, sizeof( float* ) );
	m1 = ( float** )calloc( K+1, sizeof( float* ) );
	m2 = ( float** )calloc( K+1, sizeof( float* ) );
	for( k = 0; k < K+1; k++ ){
		a[k] = ( float* )calloc( W, sizeof( float ) );
		gradient[k] = ( float* )calloc( W, sizeof( float ) );
		m1[k] = ( float* )calloc( W, sizeof( float ) );
		m2[k] = ( float* )calloc( W, sizeof( float ) );
	}
	float m1_i = 0.0f;
	float m2_i = 0.0f;

	for( k = 0; k < K+1; k++ ){

		for( j = k; j < W; j++ ){
			// get a = e^alpha
			a[k][j] = ( float )exp( alpha_[k][j] );

			// get gradients w.r.t. stochastic objective at timestep t
			gradient[k][j] = alpha_[k][j] * calcGradLogPostAlphas( alpha_[k][j], k, j );

			// update biased first moment estimate
			m1_i = beta1 * m1_i + ( 1 - beta1 ) * gradient[k][j];

			// update biased second raw moment estimate
			m2_i = beta2 * m2_i + ( 1 - beta2 ) * gradient[k][j] * gradient[k][j];

			// compute bias-corrected first moment estimate
			m1[k][j] = m1_i / ( 1 - beta1 );

			// compute bias-corrected second raw moment estimate
			m2[k][j] = m2_i / ( 1 - beta2 );

			// update parameter a due to alphas
			a[k][j] = a[k][j] - eta * m1[k][j] / ( ( float )sqrt( m2[k][j] ) + epsilon );

			// get updated alpha
			alpha_[k][j] = ( float )log( a[k][j] );
		}
	}

	// free memory
	for( k = 0; k < K+1; k++ ){
		free( gradient[k]);
		free( m1[k] );
		free( m2[k] );
	}
	free( gradient );
	free( m1 );
	free( m2 );
*/

	int k, j;
	float eta = 0.0005f;									// learning rate for a
	float** a;											// a = e^alpha
	float** gradient;									// gradient of log posterior of alpha

	a = ( float** )calloc( K+1, sizeof( float* ) );
	gradient = ( float** )calloc( K+1, sizeof( float* ) );
	for( k = 0; k < K+1; k++ ){
		a[k] = ( float* )calloc( W, sizeof( float ) );
		gradient[k] = ( float* )calloc( W, sizeof( float ) );
	}

	for( k = 0; k < K+1; k++ ){

		for( j = 0; j < W; j++ ){

			// update parameter a using a fixed learning rate
			alpha_[k][j] = alpha_[k][j] - eta * calcGradLogPostAlphas( alpha_, k, j ) ;

		}

	}

	// free memory
	for( k = 0; k < K+1; k++ ){
		free( a[k] );
		free( gradient[k]);
	}
	free( a );
	free( gradient );
}

void ModelLearning::testAlphaUpdate( float** alphas, int K, int W ){

	int k, j;

	for( k = 0; k < K+1; k++ ){

		std::cout << "k = " << k << std::endl;

		float stepsize = 0.0001f;

		for( j = 0; j < W; j++ ){

			// calculate the delta of log posterior
			alphas[k][j] += stepsize;
			float diff = calcLogPostAlphas( alphas, k );

			alphas[k][j] -= 2.0f * stepsize;
			diff -= calcLogPostAlphas( alphas, k );
			diff /= ( 2.0f * stepsize );

			// recover alphas
			alphas[k][j] += stepsize;

			// print out the gradient
			std::cout << "alpha = " << alphas[k][j] << '\t';
			std::cout << "analytical gradient: " << calcGradLogPostAlphas( alphas, k, j ) << '\t';
			std::cout << "numerical gradient: " << diff << std::endl;

		}

	}
	std::cout << std::endl;
}

float ModelLearning::calcLogPostAlphas( float** alpha, int k ){

	// calculate log posterior of alphas due to equation 46
	float logPosterior = 0.0f;
	float*** v = motif_->getV();
	float** v_bg = bg_->getV();
	int W = motif_->getW();
	int j, y, y2;

	for( j = 0; j < W; j++ ){

		// the first term of equation 46
		logPosterior -= 2.0f * logf( alpha[k][j] );

		// the second term of equation 46
		logPosterior -= Global::modelBeta * powf( Global::modelGamma, ( float )k ) / alpha[k][j];

		// the third term of equation 46
//		logPosterior += ( float )ipow( Y_[1], k ) * lgammaf( alpha[k][j] );

		/*
		// the forth term of equation 46
		for( y = 0; y < Y_[k+1]; y++ ){
			y2 = y % Y_[k];
			if( k == 0 ){
				// the first part
				logPosterior += lgamma( ( float )n_z_[k][y][j] + alpha[k][j] * v_bg[k][y] );


				// the second part
				logPosterior -= lgamma( alphas[j] * v_bg[k][y] );

				// the third part
				float* sumN = ( float* )calloc( W, sizeof( float ) );
				for( j = 0; j < W; j++ ){
					for( y = 0; y < Y_[1]; y++ ){
						sumN[j] += ( float )n_z_[0][y][j];
					}
				}
				logPosterior -= lgamma( sumN[j] + alphas[j] );
				free( sumN );


			} else {
				// the first part
				logPosterior += lgamma( ( float )n_z_[k][y][j] + alpha[k][j] * v[k-1][y2][j] );


				// the second part
				logPosterior -= lgamma( alphas[j] * v[k-1][y2][j] );

				// the third part
				logPosterior -= lgamma( ( float )n_z_[k][y][j-1] + alphas[j] );

			}

		}
*/
	}
	return logPosterior;
}

float ModelLearning::calcGradLogPostAlphas( float** alpha, int k, int j ){

	// calculate gradient of the log posterior of alphas due to equation 48
	float gradient = 0.0f;
	float*** v = motif_->getV();
	float** v_bg = bg_->getV();
	int y, y2;

	// the first term of equation 47/48
	gradient -= 2.0f / alpha[k][j];

	// the second term of equation 47/48
	gradient += Global::modelBeta * powf( Global::modelGamma, ( float )k ) / powf( alpha[k][j], 2.0f );

/*
	// the third term of equation 48
	for( y = 0; y < Y_[k+1]; y++ ){
		y2 = y % Y_[k];
		if( k == 0 ){
			float* sumN = ( float* )calloc( W, sizeof( float ) );
			for( j = 0; j < W; j++ ){
				for( y = 0; y < Y_[1]; y++ ){
					sumN[j] += n_z_[0][y][j];
				}
			}
			gradient += v_bg[k][y] * ( digammaf( ( float )n_z_[k][y][j] + alpha * v_bg[k][y] )
				- digammaf( alpha * v_bg[k][y] ) - digammaf( ( float )sumN[j] + alpha ) + digammaf( alpha ) );
			free( sumN );
		} else {
			gradient += v[k-1][y2][j] * ( digammaf( ( float )n_z_[k][y][j] + alpha * v[k-1][y2][j] )
				- digammaf( alpha * v[k-1][y2][j] ) - digammaf( ( float )n_z_[k][y][j-1] + alpha ) + digammaf( alpha ) );
		}
	}
*/

	// the third term of equation 47
//	gradient += ( float )ipow( Y_[1], k ) * digammaf( alpha[k][j] );

	/*
	// the forth term of equation 47
	for( y = 0; y < Y_[k+1]; y++ ){
		y2 = y % Y_[k];
		if( k == 0 ){
			// the first part
			gradient += v_bg[k][y] * digammaf( ( float )n_z_[k][y][j] + alpha[k][j] * v_bg[k][y] );


			// the second part
			gradient -= v_bg[k][y] * digammaf( alpha * v_bg[k][y] );

			// the third part
			int sumC = 0;
			for( y = 0; y < Y_[1]; y++ ){
				sumC += n_z_[0][y][j];
			}
			gradient -= digammaf( ( float )sumC + alpha );


		} else {
			// the first part
			gradient += v[k-1][y2][j] * digammaf( ( float )n_z_[k][y][j] + alpha[k][j] * v[k-1][y2][j] );


			// the second part
			gradient -= digammaf( alpha * v[k-1][y2][j] );

			// the third part
			gradient -= digammaf( ( float )n_z_[k][y][j-1] + alpha ) ;


		}
	}
	*/

	return gradient;
}

void ModelLearning::print(){

}

void ModelLearning::write(){

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
						+ Global::posSequenceBasename;

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
			ofile_r << std::scientific << r_[n][0] << ' ';		// print out the responsibility of not having a motif on the sequence
			ofile_pos << n+1 << '\t';							// print out the sequence number
			for( i = LW1; i > 0; i-- ){
				ofile_r << std::scientific << r_[n][i] << ' ';
				if( r_[n][i] >= cutoff ){
					ofile_pos << LW1-i << '\t';
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
				ofile_r << std::scientific << r_[n][i] << ' ';
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
