//
// Created by wanwan on 16.08.17.
//
#include "GibbsSampling.h"
#include "../Global/Global.h"

#include <boost/math/special_functions.hpp>			/* gamma and digamma function */
#include <boost/math/distributions/beta.hpp>		/* beta distribution */

GibbsSampling::GibbsSampling( Motif* motif, BackgroundModel* bg,
                              std::vector<Sequence*> seqs,
                              bool samplingQ,
                              float beta, float gamma,
                              bool initializeZ, bool samplingZ,
                              bool optimizeA, bool GibbsMHalphas, bool dissampleAlphas, bool verbose ){

    motif_              = motif;
    bg_                 = bg;
    q_                  = motif->getQ();
    seqs_               = seqs;
    beta_               = beta;
    gamma_              = gamma;
    initializeZ_        = initializeZ;
    samplingZ_          = samplingZ;
    samplingQ_          = samplingQ;
    optimizeA_          = optimizeA;
    GibbsMHalphas_      = GibbsMHalphas;
    dissampleAlphas_    = dissampleAlphas;
    verbose_            = verbose;

    // get motif (hyper-)parameters from motif class
    K_ = motif_->getK();
    W_ = motif_->getW();
    Y_ = motif_->getY();
    s_ = motif_->getS();
    A_ = motif_->getA();
    K_bg_ = ( bg_->getOrder() < K_ ) ?  bg_->getOrder() : K_;

    // allocate memory for r_[n][i], pos_[n][i], z_[n]
    r_ = ( float** )calloc( seqs_.size(), sizeof( float* ) );
    pos_ = ( float** )calloc( seqs_.size(), sizeof( float* ) );
    z_ = ( size_t* )calloc( seqs_.size(), sizeof( size_t ) );
    for( size_t n = 0; n < seqs_.size(); n++ ){
        r_[n] = ( float* )calloc( seqs_[n]->getL(), sizeof( float ) );
        pos_[n] = ( float* )calloc( seqs_[n]->getL(), sizeof( float ) );
    }

    // allocate memory for n_[k][y][j] and probs_[k][y][j]
    n_ = ( float*** )calloc( K_+1, sizeof( float** ) );
    for( size_t k = 0; k < K_+1; k++ ){
        n_[k] = ( float** )calloc( Y_[k+1], sizeof( float* ) );
        for( size_t y = 0; y < Y_[k+1]; y++ ){
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

    rngx_.seed( 42 );
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
        for( size_t y = 0; y < Y_[k+1]; y++ ){
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

    clock_t t0 = clock();

    size_t iteration = 0;

    // initialize z for all the sequences
    if( initializeZ_ ){
        EM model( motif_, bg_, seqs_ );
        // E-step: calculate posterior
        model.EStep();

        // extract initial z from the indices of the biggest responsibilities
        for( size_t n = 0; n < seqs_.size(); n++ ){
            size_t L = seqs_[n]->getL();
            float maxR = 0.f;
            size_t maxIdx = 0;
            for( size_t i = 0; i < L-W_+1; i++ ){
                if( model.getR()[n][L-W_-i] > maxR ){
                    maxR = model.getR()[n][L-W_-i];
                    maxIdx = i+1;
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
        for( size_t y = 0; y < Y_[k+1]; y++ ){
            for( size_t j = 0; j < W_; j++ ){
                n_[k][y][j] = 0.0f;
            }
        }
    }

    // 2. count k-mers for the highest order K
    for( size_t n = 0; n < seqs_.size(); n++ ){
        if( z_[n] > 0 ){
            size_t* kmer = seqs_[n]->getKmer();
            for( size_t j = 0; j < W_; j++ ){
                size_t y = kmer[z_[n]-1+j] % Y_[K_+1];
                n_[K_][y][j]++;
            }
        }
    }

    // compute k-mer counts for all the lower orders
    for( size_t k = K_; k > 0; k-- ){
        for( size_t y = 0; y < Y_[k+1]; y++ ){
            size_t y2 = y % Y_[k];
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

    // iterate over
    while( iteration < maxCGSIterations_ ){

        iteration++;

        // Collapsed Gibbs sampling position z
        if( samplingZ_ )    Collapsed_Gibbs_sampling_z();

        // Gibbs sampling fraction q
        if( samplingQ_ )	Gibbs_sample_q();

        // update alphas by stochastic optimization
        if( optimizeA_ ){

            Optimize_alphas_by_SGD_ADAM( K_, W_, eta_, iteration );

        } else if( GibbsMHalphas_ ){

            GibbsMH_sample_alphas( iteration );

            if( iteration > maxCGSIterations_ - 10 ){
                for( size_t k = 0; k < K_+1; k++ ){
                    for( size_t j = 0; j < W_; j++ ){
                        alpha_avg[k][j] += A_[k][j];
                    }
                }
            }

        } else if( dissampleAlphas_ ){

            Discrete_sample_alphas( iteration );

            if( iteration > maxCGSIterations_ - 10 ){
                for( size_t k = 0; k < K_+1; k++ ){
                    for( size_t j = 0; j < W_; j++ ){
                        alpha_avg[k][j] += A_[k][j];
                    }
                }
            }
        } else {
            std::cout << "Alphas are not optimized." << std::endl;
        }

/*
        // for making a movie out of all iterations
        // calculate probabilities
        motif_->calculateP();
        motif_->write( Global::outputDirectory,
                       Global::outputFileBasename + "_iter_" + std::to_string( iteration ) );
*/

    }

    // obtaining a motif model:
    if( GibbsMHalphas_ || dissampleAlphas_ ){
        // average alphas over the last few steps for GibbsMH
        for( size_t k = 0; k < K_+1; k++ ){
            for( size_t j = 0; j < W_; j++ ){
                A_[k][j] = alpha_avg[k][j] / 10.0f;
            }
        }
    }

    // update model parameter v
    motif_->updateV( n_, A_, K_ );

/*
    // run five steps of EM to optimize the final model with
    // the optimum model parameters v's and the fixed alphas
    EM model( motif_, bg_, seqs_, q_ );

    for( size_t step = 0; step < 5; step++ ){

        // E-step: calculate posterior
        model.EStep();

        // M-step: update model parameters
        model.MStep();
    }
*/

    // print out the optimized q
    if( verbose_ and  samplingQ_ )   std::cout << "The sampled q=" << q_ << std::endl;

    // calculate probabilities
    motif_->calculateP();

    fprintf( stdout, "\n--- Runtime for Gibbs sampling: %.4f seconds ---\n",
             ( ( float )( clock() - t0 ) ) / CLOCKS_PER_SEC );
}

void GibbsSampling::Collapsed_Gibbs_sampling_z(){

    N0_ = 0;		// reset N0

    llikelihood_ = 0.0f;

    // updated model parameters v excluding the n'th sequence
    motif_->updateV( n_, A_, K_ );

    float*** v = motif_->getV();
    float** v_bg = bg_->getV();

    // compute log odd scores s[y][j], log likelihoods of the highest order K
    motif_->calculateLinearS( v_bg, K_bg_ );

    // sampling z:
    bool remove_kmer_slowly = false;	// a flag to switch between slow and fast versions for counting k-mers

    // loop over all sequences and drop one sequence each time and update r
    for( size_t n = 0; n < seqs_.size(); n++ ){

        size_t  L = seqs_[n]->getL();
        size_t  LW1 = L - W_ + 1;
        size_t* kmer = seqs_[n]->getKmer();

        // count k-mers at position z_i+j except the n'th sequence
        // remove the k-mer counts from the sequence with the current z
        // and re-calculate the log odds scores in linear space
        /**
         * -------------- faster version of removing k-mer ------------------
         */
        float sumN = 0.0f;
        for( size_t a = 0; a < Y_[1]; a++ ){
            sumN += n_[0][a][0];
        }

        if( !remove_kmer_slowly && z_[n] > 0 ){

            for( size_t j = 0; j < W_; j++ ){

                // for k = 0:
                size_t y = kmer[z_[n]-1+j] % Y_[1];
                size_t y_bg = y % Y_[K_bg_+1];
                n_[0][y][j]--;

                v[0][y][j]= ( n_[0][y][j] + A_[0][j] * v_bg[0][y] ) / ( sumN + A_[0][j] );

                s_[y][j] = v[K_][y][j] / v_bg[K_bg_][y_bg];

                // for 1 <= k <= K_:
                for( size_t k = 1; k < K_+1; k++ ){
                    y = kmer[z_[n]-1+j] % Y_[k+1];
                    size_t y2 = y % Y_[k];
                    size_t yk = y / Y_[1];
                    y_bg = y % Y_[K_bg_+1];
                    n_[k][y][j]--;

                    if( j < K_ ){
                        v[k][y][j] = v[k-1][y2][j];
                    } else {
                        v[k][y][j] = ( n_[k][y][j] + A_[k][j] * v[k-1][y2][j] )
                                     / ( n_[k-1][yk][j-1] + A_[k][j] );
                    }

                    s_[y][j] = v[K_][y][j] / v_bg[K_bg_][y_bg];
                }

            }
        }

        /**
         * -------------- slower version of removing k-mer ------------------
         */
        if( remove_kmer_slowly ){

            // remove the k-mer counts from the sequence with the current z
            if( z_[n] > 0 ){
                for( size_t j = 0; j < W_; j++ ){
                    for( size_t k = 0; k < K_+1; k++ ){
                        size_t y = kmer[z_[n]-1+j] % Y_[k+1];
                        n_[k][y][j]--;
                    }
                }
            }

            // updated model parameters v excluding the n'th sequence
            motif_->updateV( n_, A_, K_ );
            // compute log odd scores s[y][j] of the highest order K
            motif_->calculateLinearS( bg_->getV(), K_bg_ );

        }

        /**
         * ------- sampling equation -------
         */

        // calculate responsibilities and positional priors
        // over all LW1 positions on n'th sequence:
        float normFactor = 1.0f - q_;
        float pos_i = q_ / static_cast<float>( LW1 );
        for( size_t i = 0; i < LW1; i++ ){
            pos_[n][i] = pos_i;
            r_[n][i] = 1.0f;
        }

        // todo: could be parallelized by extracting 8 sequences at once
        // ij = i+j runs over all positions in sequence
        for( size_t ij = 0; ij < LW1; ij++ ){
            // extract (K+1)-mer y from positions (i-k,...,i)
            size_t y = kmer[ij] % Y_[K_+1];
            // j runs over all motif positions
            for( size_t j = 0; j < W_; j++ ){
                r_[n][L-W_-ij+j] *= s_[y][j];
            }
        }

        // calculate responsibilities and normalize them
        for( size_t i = 0; i < LW1; i++ ){
            r_[n][i] *= pos_[n][L-W_-i];
            normFactor += r_[n][i];
        }
        for( size_t i = LW1; i < L; i++){
            r_[n][i] = 0.f;
        }
        // calculate log likelihood of sequences
        llikelihood_ += logf( normFactor );

        // normalize responsibilities and append them to a vector of posteriors
        std::vector<float> posteriors;
        // append the posterior of not having any motif on the sequence
        posteriors.push_back( ( 1.f - q_ ) / normFactor );
        for( int i = L-W_; i>=0; i-- ){
            r_[n][i] /= normFactor;
            posteriors.push_back( r_[n][i] );
        }

        // draw a new position z from the discrete distribution of posterior
        std::discrete_distribution<size_t> posterior_dist( posteriors.begin(), posteriors.end() );
        z_[n] = posterior_dist( rngx_ );

        if( z_[n] == 0 ){
            // count sequences which do not contain motifs.
            N0_++;

        } else {
            // add the k-mer counts from the current sequence with the updated z
            for( size_t j = 0; j < W_; j++ ){
                for( size_t k = 0; k < K_+1; k++ ){
                    size_t y = kmer[z_[n]-1+j] % Y_[k+1];
                    n_[k][y][j]++;
                }
            }
        }
    }
//    std::cout << "N0=" << N0_ << std::endl;
}

void GibbsSampling::Gibbs_sample_q(){

    // sampling the fraction of sequences which contain the motif
    boost::math::beta_distribution<float> q_beta_dist( ( float )seqs_.size() - ( float )N0_ + 1.0f, ( float )N0_ + 1.0f );
    q_ = quantile( q_beta_dist, ( float )rand() / ( float )RAND_MAX );

}

void GibbsSampling::Optimize_alphas_by_SGD_ADAM( size_t K, size_t W_, float eta, size_t iter){
    // update alphas using stochastic optimization algorithm ADAM
    // (DP Kingma & JL Ba 2015)

    double beta1 = 0.9;		// exponential decay rate for the moment estimates
    double beta2 = 0.999;	// exponential decay rate for the moment estimates
    double epsilon = 1e-8;	// cutoff
    double gradient;		// gradient of log posterior of alpha
    double m1;				// first moment vector (the mean)
    double m2;				// second moment vector (the uncentered variance)

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

    for( size_t k = 0; k < K_+1; k++ ){

        for( size_t j = 0; j < W_; j++ ){

            // draw 10 times in a row and take record of the last accepted sample
//		for( int step = 0; step < 10; step++ ){

            // Metropolis-Hasting sheme
            float a_prev = logf( A_[k][j] );

            float lprob_a_prev = calc_logCondProb_a( iter, a_prev, k, j );

            // draw a new 'a' from the distribution of N(a, 1)
            std::normal_distribution<float> norm_dist( a_prev,
                                                       1.0f / ( float )( k+1 ) );

            float a_new = norm_dist( rngx_ );

            float lprob_a_new = calc_logCondProb_a( iter, a_new, k, j );
            float accept_ratio;
            float uni_random;
            if( lprob_a_new < lprob_a_prev ){
                // calculate the acceptance ratio
                accept_ratio = expf( lprob_a_new - lprob_a_prev );

                // draw a random number uniformly between 0 and 1
                std::uniform_real_distribution<float> uniform_dist( 0.0f, 1.0f );
                uni_random = uniform_dist( rngx_ );

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

            A_[k][j] = expf( ( float )posterior_dist( rngx_ ) / 10.0f );
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
    gradient += beta_ * powf( gamma_, ( float )k ) / powf( A[k][j], 2.0f );

    // the third term of equation 47
    gradient += static_cast<float>( Y_[k] ) * boost::math::digamma( A[k][j] );

    // the forth term of equation 47
    for( size_t y = 0; y < Y_[k+1]; y++ ){

        if( k == 0 ){
            // the first part
            gradient += v_bg[k][y] * boost::math::digamma( n_[k][y][j] + A[k][j] * v_bg[k][y] );

            // the second part
            gradient -= v_bg[k][y] * boost::math::digamma( A[k][j] * v_bg[k][y] );

        } else {

            size_t y2 = y % Y_[k];

            // the first part
            gradient += v[k-1][y2][j] * boost::math::digamma( n_[k][y][j] + A[k][j] * v[k-1][y2][j] );

            // the second part
            gradient -= v[k-1][y2][j] * boost::math::digamma( A[k][j] * v[k-1][y2][j] );

        }

    }

    // the last term of equation 47
    for( size_t y = 0; y < Y_[k]; y++ ){

        if( k == 0 ){

            gradient -= boost::math::digamma( N + A[k][j] );

        } else if( j == 0 ){

            float sum = 0.0f;
            for( size_t a = 0; a < Y_[1]; a++ ){
                size_t ya = y * Y_[1] + a;
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
    logCondProbA -= beta_ * powf( gamma_, ( float )k ) / alpha;

    if( k == 0 ){

        for( size_t y = 0; y < Y_[k]; y++ ){

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
        for( size_t y = 0; y < Y_[k]; y++ ){

            // the third term of equation 50
            logCondProbA += boost::math::lgamma( alpha );

            // the forth term of equation 50
            for( size_t A = 0; A < Y_[1]; A++ ){

                size_t ya = y * Y_[1] + A;

                size_t y2 = ya % Y_[k];

                if( n_[k][ya][j] > 0.0f ){

                    // the first part of the forth term
                    logCondProbA += boost::math::lgamma( n_[k][ya][j] + alpha * v[k-1][y2][j] );

                    // the second part of the forth term
                    logCondProbA -= boost::math::lgamma( alpha * v[k-1][y2][j] );
                }

            }

            // the fifth term
            // Note: here it might be problematic when j = 0
            logCondProbA -= boost::math::lgamma( n_[k-1][y][j-1] + alpha );

            // !!! important: correction for the occasions when zero
            // or one k-mer is present
            if( n_[k-1][y][j-1] - 1.0f <= 0.0000001f ){

                logCondProbA = -a - beta_ * powf( gamma_, ( float )k ) / alpha;

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
        logPrior -= beta_ * powf( gamma_, ( float )k ) / alpha[k][j];

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
        logPosterior -= beta_ * powf( gamma_, ( float )k ) / A[k][j];

        //
        if( k == 0 ){

            for( size_t y = 0; y < Y_[k]; y++ ){

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
            for( size_t y = 0; y < Y_[k]; y++ ){

                // !!! important: correction for the occasions
                // when zero or one k-mer is present
                // Note: here it might be problematic when j = 0
                if( n_[k-1][y][j-1] <= 1.0f ){

                    logPosterior = - 2.0f * logf( A[k][j] ) - beta_ * powf( gamma_, ( float )k ) / A[k][j];

                } else {
                    // the third term of equation 50
                    logPosterior += boost::math::lgamma( A[k][j] );

                    // the forth term of equation 50
                    for( size_t a = 0; a < Y_[1]; a++ ){

                        size_t ya = y * Y_[1] + a;

                        size_t y2 = ya % Y_[k];

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

            for( size_t y = 0; y < Y_[k]; y++ ){

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
            for( size_t y = 0; y < Y_[k]; y++ ){

                // the third term of equation 50
                logLikelihood += boost::math::lgamma( A[k][j] );

                // the forth term of equation 50
                for( size_t a = 0; a < Y_[1]; a++ ){

                    size_t ya = y * Y_[1] + a;

                    size_t y2 = ya % Y_[k];

                    if( n_[k][ya][j] > 0.0000001f ){

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
void GibbsSampling::print(){

    // print out motif parameter v
    for( size_t j = 0; j < W_; j++ ){
        for( size_t k = 0; k < K_+1; k++ ){
            for( size_t y = 0; y < Y_[k+1]; y++ ){
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
            for( size_t y = 0; y < Y_[k+1]; y++ ){
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
