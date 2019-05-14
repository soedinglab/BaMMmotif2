//
// Created by wanwan on 16.08.17.
//
#include "GibbsSampling.h"
#include "../Global/Global.h"

#include <boost/math/special_functions.hpp>			/* gamma and digamma function */
#include <boost/math/distributions/beta.hpp>		/* beta distribution */

GibbsSampling::GibbsSampling( Motif* motif,
                              BackgroundModel* bg,
                              std::vector<Sequence*> seqs)
{

    motif_  = motif;
    bg_     = bg;
    seqs_   = seqs;
    q_      = motif->getQ();

    // get motif (hyper-)parameters from motif class
    K_      = motif_->getK();
    W_      = motif_->getW();
    s_      = motif_->getS();
    A_      = motif_->getA();
    v_      = motif_->getV();

    K_bg_   = bg_->getOrder();
    v_bg_   = bg_->getV();

    padding_= 1;

    // allocate memory for r_[n][i], pos_[n][i], z_[n]
    size_t N = seqs_.size();
    r_ = ( float** )calloc( N, sizeof( float* ) );
    pos_ = ( float** )calloc( N, sizeof( float* ) );
    z_ = ( size_t* )calloc( N, sizeof( size_t ) );
    posteriorCum_.resize( N );
    for( size_t n = 0; n < N; n++ ){
        size_t L = seqs_[n]->getL();
        r_[n] = ( float* )calloc( L+padding_+W_+1, sizeof( float ) );
        pos_[n] = ( float* )calloc( L, sizeof( float ) );
        posteriorCum_[n].resize( L-W_+2, sizeof( float ) );
    }

    // allocate memory for n_[k][y][j] and probs_[k][y][j]
    n_ = ( float*** )calloc( K_+1, sizeof( float** ) );
    for( size_t k = 0; k < K_+1; k++ ){
        n_[k] = ( float** )calloc( Global::A2powerK[k+1], sizeof( float* ) );
        for( size_t y = 0; y < Global::A2powerK[k+1]; y++ ){
            n_[k][y] = ( float* )calloc( W_, sizeof( float ) );
        }
    }

    // allocate memory for alpha_[k][j]
    m1_t_ = ( double** )calloc( K_+1, sizeof( double* ) );
    m2_t_ = ( double** )calloc( K_+1, sizeof( double* ) );
    for( size_t k = 0; k < K_+1; k++ ){
        m1_t_[k] = ( double* )calloc( W_, sizeof( double ) );
        m2_t_[k] = ( double* )calloc( W_, sizeof( double ) );
    }
    srand( 42 );
}

GibbsSampling::~GibbsSampling(){

    for( size_t n = 0; n < seqs_.size(); n++ ){
        free( r_[n] );
        free( pos_[n] );
    }
    free( r_ );
    free( pos_ );
    free( z_ );

    for( size_t k = 0; k < K_+1; k++ ){
        for( size_t y = 0; y < Global::A2powerK[k+1]; y++ ){
            free( n_[k][y] );
        }
        free( n_[k] );
        free( m1_t_[k] );
        free( m2_t_[k] );
    }
    free( n_ );
    free( m1_t_ );
    free( m2_t_ );
}


void GibbsSampling::optimize(){

    if( Global::verbose ) {
        std::cout << " __________________"  << std::endl
                  << "|*                *|" << std::endl
                  << "|  Gibbs Sampling  |" << std::endl
                  << "|*________________*|" << std::endl
                  << std::endl;
    }

    auto t0_wall = std::chrono::high_resolution_clock::now();
    // initialize z for all the sequences
    if( !Global::noInitialZ ){
        EM model( motif_, bg_, seqs_ );
        // E-step: calculate posterior
        model.EStep();
        float** r = model.getR();

        // extract initial z from the indices of the biggest responsibilities
        for( size_t n = 0; n < seqs_.size(); n++ ){
            size_t  L = seqs_[n]->getL();
            size_t  LW1 = L - W_ + 1;
            float maxR = 0.f;
            size_t maxIdx = 0;
            for( size_t i = 1; i <= LW1; i++ ){
                if( r[n][L+padding_-i] > maxR ){
                    maxR = r[n][L+padding_-i];
                    maxIdx = i;
                }
            }
            z_[n] = maxIdx;
        }

    } else {
        // initialize z with a random number
        for( size_t n = 0; n < seqs_.size(); n++ ){
            size_t LW2 = seqs_[n]->getL() - W_ + 2;
            z_[n] = static_cast<size_t>( rand() ) % LW2;
        }
    }

    // count the k-mers
    // 1. reset n_z_[k][y][j] to 0
    for( size_t k = 0; k < K_+1; k++ ){
        for( size_t y = 0; y < Global::A2powerK[k+1]; y++ ){
            for( size_t j = 0; j < W_; j++ ){
                n_[k][y][j] = 0.0f;
            }
        }
    }

    // 2. count k-mers for the highest order K
#pragma omp parallel for
    for( size_t n = 0; n < seqs_.size(); n++ ){
        if( z_[n] > 0 ){
            size_t* kmer = seqs_[n]->getKmer();
            for( size_t j = 0; j < W_; j++ ){
                size_t y = kmer[z_[n]-1+j] % Global::A2powerK[K_+1];
                n_[K_][y][j]++;
            }
        }
    }

    // compute k-mer counts for all the lower orders
    for( size_t k = K_; k > 0; k-- ){
        for( size_t y = 0; y < Global::A2powerK[k+1]; y++ ){
            size_t y2 = y % Global::A2powerK[k];
            for( size_t j = 0; j < W_; j++ ){
                n_[k-1][y2][j] += n_[k][y][j];
            }
        }
    }

/*
    // Gibbs sampling position z using empirical alphas
    for( size_t iter = 0; iter < 10; iter++ ){
        Collapsed_Gibbs_sampling_z();
    }

*/

    // vector to store the last a few alphas for certain sampling methods
    std::vector<std::vector<float>> alpha_avg( K_+1 );
    for( size_t k = 0; k < K_+1; k++ ){
        alpha_avg[k].resize( W_ );
        for( size_t j = 0; j < W_; j++ ){
            alpha_avg[k][j] = 0.0f;
        }
    }

    size_t iteration = 0;
    size_t last_steps= 10;
    // iterate over
    while( iteration < Global::maxIterations ){

        // get parameter variables with highest order before CGS
        float llikelihood_prev = llikelihood_;

        iteration++;

        // Collapsed Gibbs sampling position z
        if( !Global::noZSampling ) {
            Collapsed_Gibbs_sampling_z();
        }

        // Gibbs sampling fraction q
        if( !Global::noQSampling ){
            Gibbs_sample_q();
        }

        // update alphas by stochastic optimization
        if( !Global::noAlphaOptimization ){

            Optimize_alphas_by_SGD_ADAM( K_, W_, eta_, iteration );

        } else if( Global::GibbsMHalphas ){

            GibbsMH_sample_alphas( iteration );

            if( iteration > Global::maxIterations - last_steps ){
                for( size_t k = 0; k < K_+1; k++ ){
                    for( size_t j = 0; j < W_; j++ ){
                        alpha_avg[k][j] += A_[k][j];
                    }
                }
            }

        } else if( Global::dissampleAlphas ){

            Discrete_sample_alphas( iteration );

            if( iteration > Global::maxIterations - last_steps ){
                for( size_t k = 0; k < K_+1; k++ ){
                    for( size_t j = 0; j < W_; j++ ){
                        alpha_avg[k][j] += A_[k][j];
                    }
                }
            }

        } else {
            std::cout << "Alphas are not optimized." << std::endl;
        }

        // check the change of likelihood for convergence
        float llikelihood_diff = llikelihood_ - llikelihood_prev;

        if( Global::verbose ){
            std::cout << "it=" << iteration
                      << std::fixed
                      << "\tlog likelihood=" << llikelihood_
                      << "\tllh_diff=" << llikelihood_diff
                      << std::endl;
        }

        if( Global::makeMovie ) {
            // for making a movie out of all iterations
            // calculate probabilities
            motif_->calculateP();
            motif_->write(Global::outputDirectory,
                          Global::outputFileBasename + "_iter_"
                          + std::to_string(iteration));
        }

    }

    // obtaining a motif model:
    if( Global::GibbsMHalphas || Global::dissampleAlphas ){
        // average alphas over the last few (10) steps for GibbsMH
        for( size_t k = 0; k < K_+1; k++ ){
            for( size_t j = 0; j < W_; j++ ){
                A_[k][j] = alpha_avg[k][j] / (float)last_steps;
            }
        }
    }

    // todo: print out alphas for checking
    std::string opath = std::string( Global::outputDirectory )
                        + "/optimized.alphas" ;
    std::ofstream ofile( opath.c_str() );
    for( size_t j = 0; j < W_; j++ ){
        ofile << "pos " << std::to_string(j+1) << '\t';
    }
    ofile << '\n';
    for( size_t k = 0; k < K_+1; k++ ){
        for( size_t j = 0; j < W_; j++ ){
            ofile << std::setprecision(4) << A_[k][j] << '\t';
        }
        ofile << '\n';
    }

    // update model parameter v
    motif_->updateV( n_, A_ );

    // run 5-10 steps of EM to optimize the final model with
    // the optimum model parameters v's and the fixed alphas
    for( size_t step = 0; step < last_steps; step++ ){

        /**
         * E-step: calculate posterior
         */
        float llikelihood = 0.f;

        motif_->calculateLinearS( bg_->getV() );

        // calculate responsibilities r at all LW1 positions on sequence n
        // n runs over all sequences
#pragma omp parallel for reduction(+:llikelihood)
        for( size_t n = 0; n < seqs_.size(); n++ ){

            size_t 	L = seqs_[n]->getL();
            size_t 	LW1 = L - W_ + 1;
            size_t*	kmer = seqs_[n]->getKmer();
            float 	normFactor = 0.f;

            // initialize r_[n][i] with 1:
            // note here:   r_[n][0] is the responsibility when motif is absent
            //              and the indices of r_[n][i] is reverted, which means:
            //              r_[n][i] is the weight for motif being at position i on sequence n
            for( size_t i = 0; i <= LW1; i++ ){
                r_[n][i] = 1.0f;
            }

            // when p(z_n > 0), it runs over all positions in sequence in the region of [1, LW1]
            for( size_t i = 1; i <= LW1; i++ ){
                for( size_t j = 0; j < W_; j++ ){
                    // extract (K+1)-mer y from positions (ij-K,...,ij)
                    size_t y = kmer[i-1+j] % Global::A2powerK[K_+1];
                    r_[n][i] *= s_[y][j];
                }
            }

            // calculate the responsibilities and sum them up
            for( size_t i = 0; i <= LW1; i++ ){
                r_[n][i] *= pos_[n][i];
                normFactor += r_[n][i];
            }

            // normalize responsibilities in the region [0, LW1]
            for( size_t i = 0; i <= LW1; i++ ){
                r_[n][i] /= normFactor;
            }

            // calculate log likelihood over all sequences
            llikelihood += logf( normFactor );
        }

        llikelihood_ = llikelihood;

        /**
         * M-step: update model parameters
         */
        // reset the fractional counts n
        for( size_t k = 0; k < K_+1; k++ ){
            for( size_t y = 0; y < Global::A2powerK[k+1]; y++ ){
                for( size_t j = 0; j < W_; j++ ){
                    n_[k][y][j] = 0.0f;
                }
            }
        }

        // compute fractional occurrence counts for the highest order K
        // n runs over all sequences
#pragma omp parallel for
        for( size_t n = 0; n < seqs_.size(); n++ ){
            size_t  L = seqs_[n]->getL();
            size_t  LW1 = L - W_ + 1;
            size_t* kmer = seqs_[n]->getKmer();

            // it runs over all positions i on sequence n in the region of [1, LW1]
            for( size_t i = 1; i <= LW1; i++ ){
                for (size_t j = 0; j < W_; j++) {
                    size_t y = kmer[i-1+j] % Global::A2powerK[K_+1];
                    //n_[K_][y][j] += r_[n][i];
                    atomic_float_add(&(n_[K_][y][j]), r_[n][i]);
                }
            }
        }

        // compute fractional occurrence counts from higher to lower order
        // k runs over all lower orders
        for( size_t k = K_; k > 0; k-- ){
            for( size_t y = 0; y < Global::A2powerK[k+1]; y++ ){
                size_t y2 = y % Global::A2powerK[k];
                for( size_t j = 0; j < W_; j++ ){
                    n_[k-1][y2][j] += n_[k][y][j];
                }
            }
        }

        // update model parameters v[k][y][j] with updated k-mer counts, alphas and model order
        motif_->updateV( n_, A_ );
    }


    // print out the optimized q
    if( Global::verbose and !Global::noQSampling ){
        std::cout << "The sampled q=" << q_ << std::endl;
    }

    // calculate probabilities
    motif_->calculateP();

    auto t1_wall = std::chrono::high_resolution_clock::now();
    auto t_diff = std::chrono::duration_cast<std::chrono::duration<double>>(t1_wall-t0_wall);
    std::cout << "\n--- Runtime for Gibbs sampling: " << t_diff.count() << " seconds ---\n";

}

void GibbsSampling::Collapsed_Gibbs_sampling_z(){

    N0_ = 0;		// reset N0

    llikelihood_ = 0.0f;

    // updated model parameters v excluding the n'th sequence
    motif_->updateV( n_, A_ );

    // compute log odd scores s[y][j], log likelihoods of the highest order K
    motif_->calculateLinearS( v_bg_ );

    // reset k_bg
    // todo: check if this is correct
    size_t k_bg = ( K_ > K_bg_ ) ? K_bg_ : K_;
    // sampling z:
//    bool remove_kmer_slowly = false;	// a flag to switch between slow and fast versions for counting k-mers

    // loop over all sequences and drop one sequence each time and update r
//#pragma omp parallel for
    // todo: this parallelization does not give the same result yet
    for( size_t n = 0; n < seqs_.size(); n++ ){

        size_t  L   = seqs_[n]->getL();
        size_t  LW1 = L - W_ + 1;
        size_t* kmer = seqs_[n]->getKmer();

        // count k-mers at position z_i+j except the n'th sequence
        // remove the k-mer counts from the sequence with the current z
        // and re-calculate the log odds scores in linear space
        /**
         * -------------- faster version of removing k-mer ------------------
         */
        float sumN = 0.0f;
        for( size_t a = 0; a < Global::A2powerK[1]; a++ ){
            sumN += n_[0][a][0];
        }

        if( /* !remove_kmer_slowly && */ z_[n] > 0 ){

            for( size_t j = 0; j < W_; j++ ){

                // for k = 0:
                size_t y = kmer[z_[n]-1+j] % Global::A2powerK[1];
                size_t y_bg = y % Global::A2powerK[k_bg+1];
                n_[0][y][j]--;

                v_[0][y][j]= ( n_[0][y][j] + A_[0][j] * v_bg_[0][y] )
                            / ( sumN + A_[0][j] );

                s_[y][j] = v_[K_][y][j] / v_bg_[k_bg][y_bg];

                // for 1 <= k <= K_:
                for( size_t k = 1; k < K_+1; k++ ){
                    y = kmer[z_[n]-1+j] % Global::A2powerK[k+1];
                    size_t y2 = y % Global::A2powerK[k];
                    size_t yk = y / Global::A2powerK[1];
                    y_bg = y % Global::A2powerK[k_bg+1];
                    n_[k][y][j]--;

                    if( j < K_ ){
                        v_[k][y][j] = v_[k-1][y2][j];
                    } else {
                        v_[k][y][j] = ( n_[k][y][j] + A_[k][j] * v_[k-1][y2][j] )
                                     / ( n_[k-1][yk][j-1] + A_[k][j] );
                    }

                    s_[y][j] = v_[K_][y][j] / v_bg_[k_bg][y_bg];
                }
            }
        }

        /**
         * -------------- slower version of removing k-mer ------------------
         */
/*
        if( remove_kmer_slowly ){

            // remove the k-mer counts from the sequence with the current z
            if( z_[n] > 0 ){
                for( size_t j = 0; j < W_; j++ ){
                    for( size_t k = 0; k < K_+1; k++ ){
                        size_t y = kmer[z_[n]-1+j] % Global::A2powerK[k+1];
                        n_[k][y][j]--;
                    }
                }
            }

            // updated model parameters v excluding the n'th sequence
            motif_->updateV( n_, A_ );
            // compute log odd scores s[y][j] of the highest order K
            motif_->calculateLinearS( bg_->getV() );
        }
*/

        /**
         * ------- sampling equation -------
         */
        // calculate responsibilities and positional priors
        // over all LW1 positions on n'th sequence:
        float normFactor = 0.f;
        float pos_i = q_ / static_cast<float>( LW1 );
        for( size_t i = 1; i <= LW1; i++ ){
            pos_[n][i] = pos_i;
            r_[n][L+padding_-i] = 1.0f;
        }

        // todo: could be parallelized by extracting 8 sequences at once
        // ij = i+j runs over all positions in sequence
        for( size_t ij = 0; ij < L; ij++ ){
            // extract (K+1)-mer y from positions (i-k,...,i)
            size_t y = kmer[ij] % Global::A2powerK[K_+1];
            // j runs over all motif positions
            for( size_t j = 0; j < W_; j++ ){
                // todo: time-consuming!
                r_[n][L+padding_-ij+j-1] *= s_[y][j];
            }
        }

        // calculate responsibilities and normalize them
        posteriorCum_[n][0] = r_[n][0] = pos_[n][0];

        normFactor += r_[n][0];
        for( size_t i = 1; i <= LW1; i++ ){
            r_[n][L+padding_-i] *= pos_[n][i];
            posteriorCum_[n][i] = r_[n][L+padding_-i]+posteriorCum_[n][i-1];
            normFactor += r_[n][L+padding_-i];
        }

        // calculate log likelihood of sequences
        llikelihood_ += logf( normFactor );


/*
        // normalize responsibilities and append them to a vector of posteriors
        std::vector<float> posteriors;
        // append the posterior of not having any motif on the sequence
        posteriors.push_back( r_[n][0] );
        for( size_t i = 1; i <= LW1; i++ ){
            // todo: slightly time consuming!
            posteriors.push_back( r_[n][L+padding_-i] );
        }

        // draw a new position z from the discrete distribution of posterior
        // todo: very time-consuming!
        std::discrete_distribution<size_t> posterior_dist( posteriors.begin(),
                                                           posteriors.end() );

        z_[n] = posterior_dist( Global::rngx );
*/

        std::uniform_real_distribution<float> uni_dist(0.f, normFactor);
        z_[n] = (size_t) std::distance( posteriorCum_[n].begin(),
                                        std::lower_bound( posteriorCum_[n].begin(),
                                                          posteriorCum_[n].end(),
                                                          uni_dist( Global::rngx ) ));
        if( z_[n] == 0 ){
            // count sequences which do not contain motifs.
            N0_++;

        } else {
            // add the k-mer counts from the current sequence with the updated z
            for( size_t j = 0; j < W_; j++ ){
                for( size_t k = 0; k < K_+1; k++ ){
                    size_t y = kmer[z_[n]-1+j] % Global::A2powerK[k+1];
                    n_[k][y][j]++;
                }
            }
        }
    }

    if(Global::verbose){
        std::cout << "N0=" << N0_ << std::endl;
    }
}

void GibbsSampling::Gibbs_sample_q(){

    // sampling the fraction of sequences which contain the motif
    boost::math::beta_distribution<float> q_beta_dist( ( float )seqs_.size() - ( float )N0_ + 1.0f, ( float )N0_ + 1.0f );
    q_ = quantile( q_beta_dist, ( float )rand() / ( float )RAND_MAX );

}

void GibbsSampling::Optimize_alphas_by_SGD_ADAM( size_t K,
                                                 size_t W_,
                                                 float eta,
                                                 size_t iter)
{
    // update alphas using stochastic optimization algorithm ADAM
    // (DP Kingma & JL Ba 2015)

    double beta1 = 0.9;		// exponential decay rate for the moment estimates
    double beta2 = 0.999;	// exponential decay rate for the moment estimates
    double epsilon = 1e-8;	// cutoff
    double gradient;		// gradient of log posterior of alpha
    double m1;				// first moment vector (the mean)
    double m2;				// second moment vector (the un-centered variance)

    double t = static_cast<double>( iter );

    for( size_t k = 0; k < K+1; k++ ){

        for( size_t j = 0; j < W_; j++ ){

            // re-parameterise alpha on log scale: alpha = e^a
            float a = logf( A_[k][j] );

            // get gradients w.r.t. stochastic objective at timestep t
            gradient = A_[k][j] * calc_gradient_alphas( A_, k, j );

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
            a += static_cast<float>( eta * m1 / ( ( sqrt( m2 ) + epsilon ) * sqrt( t ) ) );

            A_[k][j] = expf( a );
        }

    }

}

void GibbsSampling::GibbsMH_sample_alphas( size_t iter ){
    // Gibbs sampling alphas in exponential space
    // with Metropolis-Hastings algorithm
    // build a uniform distribution between 0 and 1
    std::uniform_real_distribution<float> uniform_dist( 0.0f, 1.0f );

    for( size_t k = 0; k < K_+1; k++ ){

        for( size_t j = 0; j < W_; j++ ){

            // draw 10 times in a row and take record of the last accepted sample
//		for( int step = 0; step < 10; step++ ){

            // Metropolis-Hasting scheme
            float a_prev = logf( A_[k][j] );

            float lprob_a_prev = calc_logCondProb_a( iter, a_prev, k, j );

            // draw a new 'a' from the distribution of N(a, 1)
            std::normal_distribution<float> norm_dist( a_prev,
                                                       1.0f / ( float )( k+1 ) );

            float a_new = norm_dist( Global::rngx );

            float lprob_a_new = calc_logCondProb_a( iter, a_new, k, j );
            float accept_ratio;
            float uni_random;
            if( lprob_a_new < lprob_a_prev ){
                // calculate the acceptance ratio
                accept_ratio = expf( lprob_a_new - lprob_a_prev );

                // draw a random number uniformly between 0 and 1
                uni_random = uniform_dist( Global::rngx );

                // accept the trial sample if the ratio is not smaller than
                // a random number between (0,1)
                if( accept_ratio >= uni_random ){
                    A_[k][j] = expf( a_new );
                }

            } else {
                // accept the trial sample
                A_[k][j] = expf( a_new );

            }
        }
//		}
    }
}

void GibbsSampling::Discrete_sample_alphas( size_t iter ){
    // sample an alpha from the discrete distribution of its log posterior

    for( size_t k = 0; k < K_+1; k++ ){

        for( size_t j = 0; j < W_; j++ ){

            std::vector<float> condProb;

            float condProb_new;

            float base = calc_logCondProb_a( iter, 0.0, k, j );

            for( size_t it = 0; it < 100; it++ ){

                condProb_new = expf( calc_logCondProb_a( iter, ( float )it / 10.0f, k, j ) - base );

                condProb.push_back( condProb_new );

            }

            std::discrete_distribution<> posterior_dist( condProb.begin(),
                                                         condProb.end() );

            A_[k][j] = expf( ( float )posterior_dist( Global::rngx ) / 10.0f );
        }
    }
}

float GibbsSampling::calc_gradient_alphas( float** A, size_t k, size_t j ){
    // calculate partial gradient of the log posterior of alphas
    // due to equation 47 in the theory
    // Note that j >= k

    float 		gradient = 0.0f;
    float*** 	v = motif_->getV();
    float** 	v_bg = bg_->getV();
    float 		N = static_cast<float>( seqs_.size() - 1 );

    // the first term of equation 47
    gradient -= 2.0f / A[k][j];

    // the second term of equation 47
    gradient += Global::modelBeta * powf( Global::modelGamma, ( float )k ) / powf( A[k][j], 2.0f );

    // the third term of equation 47
    gradient += static_cast<float>( Global::A2powerK[k] ) * boost::math::digamma( A[k][j] );

    // the forth term of equation 47
    for( size_t y = 0; y < Global::A2powerK[k+1]; y++ ){

        if( k == 0 ){
            // the first part
            gradient += v_bg[k][y] * boost::math::digamma( n_[k][y][j] + A[k][j] * v_bg[k][y] );

            // the second part
            gradient -= v_bg[k][y] * boost::math::digamma( A[k][j] * v_bg[k][y] );

        } else {

            size_t y2 = y % Global::A2powerK[k];

            // the first part
            gradient += v[k-1][y2][j] * boost::math::digamma( n_[k][y][j] + A[k][j] * v[k-1][y2][j] );

            // the second part
            gradient -= v[k-1][y2][j] * boost::math::digamma( A[k][j] * v[k-1][y2][j] );

        }

    }

    // the last term of equation 47
    for( size_t y = 0; y < Global::A2powerK[k]; y++ ){

        if( k == 0 ){

            gradient -= boost::math::digamma( N + A[k][j] );

        } else if( j == 0 ){

            float sum = 0.0f;
            for( size_t a = 0; a < Global::A2powerK[1]; a++ ){
                size_t ya = y * Global::A2powerK[1] + a;
                sum += n_[k][ya][j];
            }
            gradient -= boost::math::digamma( sum + A[k][j] );

        } else {

            gradient -= boost::math::digamma( n_[k-1][y][j-1] + A[k][j] );

        }
    }
    return gradient;
}

float GibbsSampling::calc_logCondProb_a( size_t iteration, float a, size_t k, size_t j ){
    // calculate partial log conditional probabilities of a's
    // due to equation 50 in the theory

    float logCondProbA = 0.0f;
    float*** v = motif_->getV();
    float** v_bg = bg_->getV();
    float N = static_cast<float>( seqs_.size() - 1 );

    // get alpha by alpha = e^a
    float alpha = expf( a );

    // the first term of equation 50
    logCondProbA -= a;

    // the second term of equation 50
    logCondProbA -= Global::modelBeta * powf( Global::modelGamma, ( float )k ) / alpha;

    if( k == 0 or j==0 ){

        for( size_t y = 0; y < Global::A2powerK[k]; y++ ){

            // the third term of equation 50
            logCondProbA += boost::math::lgamma( alpha );

            // the first part of the forth term of equation 50
            logCondProbA += boost::math::lgamma( n_[k][y][j] + alpha * v_bg[k][y] );

            // the second part of the forth term of equation 50
            logCondProbA -= boost::math::lgamma( alpha * v_bg[k][y] );

        }

        // the fifth term of equation 50
        logCondProbA -= boost::math::lgamma( N + alpha );

    } else {

        // the third, forth and fifth terms of equation 50
        for( size_t y = 0; y < Global::A2powerK[k]; y++ ){

            // the third term of equation 50
            logCondProbA += boost::math::lgamma( alpha );

            // the forth term of equation 50
            for( size_t A = 0; A < Global::A2powerK[1]; A++ ){

                size_t ya = y * Global::A2powerK[1] + A;

                size_t y2 = ya % Global::A2powerK[k];

                if( n_[k][ya][j] > 0.0f ){

                    //todo: debug: Invalid read of size 4
                    // the first part of the forth term
                    logCondProbA += boost::math::lgamma( n_[k][ya][j] + alpha * v[k-1][y2][j] );

                    // the second part of the forth term
                    logCondProbA -= boost::math::lgamma( alpha * v[k-1][y2][j] );
                }

            }

            // the fifth term
            // Note: here it might be problematic when j = 0
            // todo: debug: problem when j = 0
            logCondProbA -= boost::math::lgamma( n_[k-1][y][j-1] + alpha );

            // !!! important: correction for the occasions when zero
            // or one k-mer is present
            if( n_[k-1][y][j-1] - 1.0f <= Global::eps ){

                logCondProbA = - a - Global::modelBeta * powf( Global::modelGamma, ( float )k ) / alpha;

            }

        }

    }

    return logCondProbA;
}

float GibbsSampling::calc_prior_alphas( float** alpha, size_t k ){
    // calculate partial log conditional probabilities of alphas
    // due to equation 34 in the theory

    float logPrior = 0.0f;
    for( size_t j = 0; j < W_; j++ ){

        // the first term of equation 46
        logPrior -= 2.0f * logf( alpha[k][j] );

        // the second term of equation 46
        logPrior -= Global::modelBeta * powf( Global::modelGamma, ( float )k )
                    / alpha[k][j];

    }

    return logPrior;
}

float GibbsSampling::calc_lposterior_alphas( float** A, size_t k ){
    // calculate partial log posterior of alphas for the order k
    // due to equation 50 in the theory

    float logPosterior = 0.0f;
    float*** v = motif_->getV();
    float** v_bg = bg_->getV();
    float N = static_cast<float>( seqs_.size() - 1 );

    for( size_t j = 0; j < W_; j++ ){

        // the prior
        logPosterior -= 2.0f * logf( A[k][j] );

        // the prior
        logPosterior -= Global::modelBeta * powf( Global::modelGamma,
                                                  ( float )k ) / A[k][j];

        //
        if( k == 0 ){

            for( size_t y = 0; y < Global::A2powerK[k]; y++ ){

                // the third term of equation 50
                logPosterior += boost::math::lgamma( A[k][j] );

                // the first part of the forth term of equation 50
                logPosterior += boost::math::lgamma( n_[k][y][j] + A[k][j] * v_bg[k][y] );

                // the second part of the forth term of equation 50
                logPosterior -= boost::math::lgamma( A[k][j] * v_bg[k][y] );

            }

            // the fifth term of equation 50
            logPosterior -= boost::math::lgamma( N + A[k][j] );

        } else {

            // the third, forth and fifth terms of equation 50
            for( size_t y = 0; y < Global::A2powerK[k]; y++ ){

                // !!! important: correction for the occasions
                // when zero or one k-mer is present
                // Note: here it might be problematic when j = 0
                if( n_[k-1][y][j-1] <= 1.0f ){

                    logPosterior = - 2.0f * logf( A[k][j] ) -
                            Global::modelBeta * powf( Global::modelGamma, ( float )k ) / A[k][j];

                } else {
                    // the third term of equation 50
                    logPosterior += boost::math::lgamma( A[k][j] );

                    // the forth term of equation 50
                    for( size_t a = 0; a < Global::A2powerK[1]; a++ ){

                        size_t ya = y * Global::A2powerK[1] + a;

                        size_t y2 = ya % Global::A2powerK[k];

                        if( n_[k][ya][j] > 0 ){

                            // the first part of the forth term
                            logPosterior += boost::math::lgamma( n_[k][ya][j] + A[k][j] * v[k-1][y2][j] );

                            // the second part of the forth term
                            logPosterior -= boost::math::lgamma( A[k][j] * v[k-1][y2][j] );
                        }

                    }

                    // the fifth term
                    logPosterior -= boost::math::lgamma( n_[k-1][y][j-1] + A[k][j] );
                }
            }
        }
    }
    return logPosterior;
}

float GibbsSampling::calc_llikelihood_alphas( float** A, size_t k ){
    // calculate partial log likelihood of alphas
    // due to equation 50 in the theory

    float logLikelihood = 0.0f;
    float*** v = motif_->getV();
    float** v_bg = bg_->getV();
    float N = static_cast<float>( seqs_.size() - 1 );

    if( k == 0 ){

        for( size_t j = 0; j < W_; j++ ){

            for( size_t y = 0; y < Global::A2powerK[k]; y++ ){

                // the third term of equation 50
                logLikelihood += boost::math::lgamma( A[k][j] );

                // the first part of the forth term of equation 50
                logLikelihood += boost::math::lgamma( n_[k][y][j] + A[k][j] * v_bg[k][y] );

                // the second part of the forth term of equation 50
                logLikelihood -= boost::math::lgamma( A[k][j] * v_bg[k][y] );

            }

            // the fifth term of equation 50
            logLikelihood -= boost::math::lgamma( N + A[k][j] );
        }

    } else {

        for( size_t j = 0; j < W_; j++ ){
            // the third, forth and fifth terms of equation 50
            for( size_t y = 0; y < Global::A2powerK[k]; y++ ){

                // the third term of equation 50
                logLikelihood += boost::math::lgamma( A[k][j] );

                // the forth term of equation 50
                for( size_t a = 0; a < Global::A2powerK[1]; a++ ){

                    size_t ya = y * Global::A2powerK[1] + a;

                    size_t y2 = ya % Global::A2powerK[k];

                    if( n_[k][ya][j] > 1.e-6f ){

                        // the first part of the forth term
                        logLikelihood += boost::math::lgamma( n_[k][ya][j] + A[k][j] * v[k-1][y2][j] );

                        // the second part of the forth term
                        logLikelihood -= boost::math::lgamma( A[k][j] * v[k-1][y2][j] );
                    }

                }

                // the fifth term
                // Note: here it might be problematic when j = 0
                logLikelihood -= boost::math::lgamma( n_[k-1][y][j-1] + A[k][j] );

            }
        }
    }

    return logLikelihood;
}

float GibbsSampling::getQ(){
    return q_;
}

void GibbsSampling::printZ(){
    for(size_t n = 0; n < seqs_.size(); n++ ){
        std::cout << z_[n] << std::endl;
    }
}

void GibbsSampling::print(){

    // print out motif parameter v
    for( size_t j = 0; j < W_; j++ ){
        for( size_t k = 0; k < K_+1; k++ ){
            for( size_t y = 0; y < Global::A2powerK[k+1]; y++ ){
                std::cout << std::setprecision(5) <<
                          motif_->getV()[k][y][j] << '\t';
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

}

void GibbsSampling::write( char* odir, std::string basename, bool ss ){

    /**
     * 	 * save BaMM (hyper-)parameters in flat files:
     * (1) posSequenceBasename.counts:	    refined fractional counts of (k+1)-mers
     * (2) posSequenceBasename.weights:     responsibilities, posterior distributions
     * (3) posSequenceBasename.alphas:	    optimized hyper-parameter alphas
     * (4) posSequenceBasename.positions:   position of motif pattern on each sequence
     */

    std::string opath = std::string( odir ) + '/' + basename;

    // output (k+1)-mer counts: n[k][y][j]
    std::string opath_n = opath + ".counts";
    std::ofstream ofile_n( opath_n.c_str() );
    for( size_t j = 0; j < W_; j++ ){
        for( size_t k = 0; k < K_+1; k++ ){
            for( size_t y = 0; y < Global::A2powerK[k+1]; y++ ){
                ofile_n << static_cast<int>( n_[k][y][j] ) << '\t';
            }
            ofile_n << std::endl;
        }
        ofile_n << std::endl;
    }

    // output position(s) of motif(s): pos_[n][i]
    std::string opath_pos = opath + ".positions";
    std::ofstream ofile_pos( opath_pos.c_str() );

    ofile_pos << "seq\tlength\tstrand\tstart..end\tmotif" << std::endl;

    for( size_t n = 0; n < seqs_.size(); n++ ){

        size_t seqlen = seqs_[n]->getL();
        if( !ss )	seqlen = ( seqlen - 1 ) / 2;

        if( z_[n] > 0 ){
            ofile_pos << seqs_[n]->getHeader() << '\t'
                      << seqlen << '\t'
                      << ( ( z_[n] < seqlen ) ? '+' : '-' ) << '\t'
                      << z_[n] << ".." << z_[n]+W_-1 << '\t';
            for( size_t b = 0; b < W_; b++ ){
                ofile_pos << Alphabet::getBase( seqs_[n]->getSequence()[z_[n]+b-1] );
            }
            ofile_pos << std::endl;
        }
    }

    // output hyper-parameter alphas: alpha[k][j]
    std::string opath_alpha = opath + ".alphas";
    std::ofstream ofile_alpha( opath_alpha.c_str() );
    for( size_t k = 0; k < K_+1; k++ ){
        ofile_alpha << "> k=" << k << std::endl;
        for( size_t j = 0; j < W_; j++ ){
            ofile_alpha << std::setprecision( 3 ) << A_[k][j] << '\t';
        }
        ofile_alpha << std::endl;
    }

    // output responsibilities r[n][i]
    bool write_r = false;		// flag for writing out the responsibilities
    if( write_r ){
        std::string opath_r = opath + ".weights";
        std::ofstream ofile_r( opath_r.c_str() );
        for( size_t n = 0; n < seqs_.size(); n++ ){
            ofile_r << std::scientific << std::setprecision( 2 )
                    << r_[n][0] << ' ';
            size_t LW1 = seqs_[n]->getL() - W_ + 1;
            for( size_t i = LW1; i > 0; i-- ){
                ofile_r << std::setprecision( 2 ) << r_[n][i] << ' ';
            }
            ofile_r << std::endl;
        }
    }
}
