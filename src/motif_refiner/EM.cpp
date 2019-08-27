//
// Created by wanwan on 16.08.17.
//
#include "EM.h"

EM::EM( Motif* motif, BackgroundModel* bgModel,
        std::vector<Sequence*> seqs ){

    motif_      = motif;
	bgModel_    = bgModel;
	seqs_       = seqs;

    // get motif (hyper-)parameters from motif class
	K_          = Global::modelOrder;
	W_          = motif_->getW();
	s_          = motif_->getS();
	A_          = motif_->getA();
    q_          = motif_->getQ();   // get fractional prior q from initial motifs

    // allocate memory for r_[n][i], prior_[n][i], z_[n]
    padding_    = 1;                // add to reserve enough memory in r for fast EM version
	r_          = ( float** )calloc( seqs_.size(), sizeof( float* ) );
	prior_      = ( float** )calloc( seqs_.size(), sizeof( float* ) );
	for( size_t n = 0; n < seqs_.size(); n++ ){
		r_[n]   = ( float* )calloc( seqs_[n]->getL()+padding_+W_+1, sizeof( float ) );
		prior_[n] = ( float* )calloc( seqs_[n]->getL()+1, sizeof( float ) );
	}

	// allocate memory for n_[k][y][j]
	n_ = ( float*** )calloc( K_+1, sizeof( float** ) );
	for( size_t k = 0; k < K_+1; k++ ){
		n_[k] = ( float** )calloc( Global::A2powerK[k+1], sizeof( float* ) );
		for( size_t y = 0; y < Global::A2powerK[k+1]; y++ ){
			n_[k][y] = ( float* )calloc( W_, sizeof( float ) );
		}
	}

    // todo: for positional prior
    beta1_      = 5;                // a hyper-parameter for estimating the positional prior
    beta2_      = 1000;             // for method 2 and 3
    N0_         = 0.f;

    size_t LW1_ = seqs_[0]->getL()-W_+1;
    pi_         = ( float* )calloc( LW1_+1, sizeof( float ) );
    b_vector_   = Eigen::VectorXf::Zero( LW1_ );
    si_         = Eigen::VectorXf::Zero( LW1_ );
    Ni_         = Eigen::VectorXf::Zero( LW1_ );
    A_matrix_   = getAmatrix( LW1_ );

}

EM::~EM(){
    for( size_t n = 0; n < seqs_.size(); n++ ){
        free( r_[n] );
        free( prior_[n] );
    }
    free( r_ );
    free( prior_ );

    for( size_t k = 0; k < K_+1; k++ ){
        for( size_t y = 0; y < Global::A2powerK[k+1]; y++ ){
            free( n_[k][y] );
        }
        free( n_[k] );
    }
    free( n_ );

    free( pi_ );
}

int EM::optimize(){

    auto    t0_wall = std::chrono::high_resolution_clock::now();

    bool 	iterate = true;

    // store model probs before each optimization step
    std::vector<std::vector<float>> v_before;
    v_before.resize( Global::A2powerK[K_+1] );
    for( size_t y = 0; y < Global::A2powerK[K_+1]; y++ ){
        v_before[y].resize( W_ );
    }

    // initialized positional priors
    updatePrior();

    // todo: for positional prior
    // todo:=====================================================
    if( Global::optimizePos ) {
        // pre-define hyper-parameters for optimizing position priors
        size_t LW1 = seqs_[0]->getL()-W_+1;

        norm_ = 0.0f;
        for (size_t i = 1; i <= LW1; i++) {
            std::normal_distribution<> nd{LW1 / 2.0f, (float) LW1};
            si_[i-1] = nd( Global::rngx ) / LW1;
            pi_[i] = expf(si_[i-1]);
            norm_ += pi_[i];
        }
        for (size_t i = 1; i <= LW1; i++) {
            pi_[i] /= norm_;
        }
    }
    // todo:=====================================================

    // iterate over
    size_t iteration = 0;
    while( iterate && ( iteration < Global::maxIterations ) ){

        // get parameter variables with highest order before EM
        float llikelihood_prev = llikelihood_;

        for( size_t y = 0; y < Global::A2powerK[K_+1]; y++ ){
            for( size_t j = 0; j < W_; j++ ){
                v_before[y][j] = motif_->getV()[K_][y][j];
            }
        }

        if( !Global::slowEM ){
            EStep();            // E-step: calculate posterior
            MStep();            // M-step: update model parameters
        } else {
            EStep_slow();
            MStep_slow();
        }

        // optimize hyper-parameter q in the first step
        if( Global::optimizeQ and iteration <= 0 ) {
            optimizeQ();
            updatePrior();
        }

        // todo: optimize positional prior pos
        // todo:=====================================================
        // note: this only works for sequences with the same length
        if( Global::optimizePos /*and iteration <= 0*/ ){
            optimizePos();
        }
        // todo:=====================================================

        // check parameter difference for convergence
        float v_diff = 0.0f;
        for( size_t y = 0; y < Global::A2powerK[K_+1]; y++ ){
            for( size_t j = 0; j < W_; j++ ){
                v_diff += fabsf( motif_->getV()[K_][y][j] - v_before[y][j] );
            }
        }

        // check the change of likelihood for convergence
        float llikelihood_diff = llikelihood_ - llikelihood_prev;

        if( Global::verbose ) {
            std::cout << "i=" << iteration
                      << std::fixed
                      << "\tllh=" << llikelihood_
                      << "\tllh_diff=" << llikelihood_diff
                      << "\t\tmodel_diff=" << v_diff
                      << std::endl;
        }

        // stopping criteria
        if( ( llikelihood_diff < Global::EMepsilon or v_diff < Global::EMepsilon )
            and iteration > 10 ){
            iterate = false;
        }

        // todo: for making a movie out of all iterations
        if( Global::makeMovie ) {
            // calculate probabilities
            motif_->calculateP();
            // write out the learned model
            std::string opath = std::string( Global::outputDirectory ) + "/movie_factory" ;
            char *odir = new char[opath.length() + 1];
            std::strcpy(odir, opath.c_str());
            createDirectory( odir );
            motif_->write(odir, "movie_clap" + std::to_string(iteration));
        }

        // todo:
        if( Global::makeMovie and Global::savePi ){
            // pre-define hyper-parameters for optimizing position priors
            size_t LW1 = seqs_[0]->getL()-W_+1;
            std::string opath = std::string( Global::outputDirectory ) +'/'
                                + Global::outputFileBasename
                                + "_movie_clap_" + std::to_string(iteration) + ".pi";
            std::ofstream ofile( opath.c_str() );;
            for (size_t i = 1; i <= LW1; i++) {
                ofile << pi_[i] << std::endl;
            }
        }

        // todo: print out for checking
        if( Global::debug ) {
            motif_->calculateP();
            printR();
        }

        iteration++;
    }

    // update probabilities
    motif_->calculateP();

    auto t1_wall = std::chrono::high_resolution_clock::now();
    auto t_diff = std::chrono::duration_cast<std::chrono::duration<double>>(t1_wall-t0_wall);
    std::cout << "\n--- Runtime for EM: " << t_diff.count()
              << " seconds (" << std::to_string(iteration)
              << " iterations) ---\n";

    return 0;
}

void EM::EStep(){

    float llikelihood = 0.0f;

    motif_->calculateLinearS( bgModel_->getV() );

    // calculate responsibilities r at all LW1 positions on sequence n
    // n runs over all sequences
#pragma omp parallel for reduction(+:llikelihood)
    for( size_t n = 0; n < seqs_.size(); n++ ){

        size_t 	L = seqs_[n]->getL();
        size_t  LW1 = L - W_ + 1;
        size_t*	kmer = seqs_[n]->getKmer();
        float 	normFactor = 0.f;

        // initialize r_[n][i] and prior_[n][i]:
        // note here:   r_[n][0] and prior_[n][0] are for motif is absent
        //              and the indices of r_[n][i] is reverted, which means:
        //              r_[n][i] is the weight for motif being at position L-i on sequence n
        //              the purpose of doing the reverting is to improve the computation speed

        for( size_t i = 1; i <= LW1; i++ ){
            r_[n][L+padding_-i] = 1.0f;
        }

        // when p(z_n > 0), ij = i+j runs over all positions in sequence
        // r index goes from (0, L]

        for( size_t ij = 0; ij < L; ij++ ){
            // extract (K+1)-mer y from positions (ij-K,...,ij)
            size_t y = kmer[ij] % Global::A2powerK[K_+1];
            for( size_t j = 0; j < W_; j++ ){
                r_[n][L+padding_-ij+j-1] *= s_[y][j];
            }
        }

        // calculate the responsibilities and sum them up
        r_[n][0] = prior_[n][0];
        normFactor += r_[n][0];
        for( size_t i = 1; i <= LW1; i++ ){
            r_[n][L+padding_-i] *= prior_[n][i];
            normFactor += r_[n][L+padding_-i];
        }

        // normalize responsibilities
        r_[n][0] /= normFactor;
        for( size_t i = 1; i <= LW1; i++ ){
            r_[n][L+padding_-i] /= normFactor;
        }

        // calculate log likelihood over all sequences
        llikelihood += logf( normFactor );
    }

    llikelihood_ = llikelihood;
}

void EM::MStep(){

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
        size_t* kmer = seqs_[n]->getKmer();

        // ij = i+j runs over all positions i on sequence n
        // r index goes from (0, L]
        for( size_t ij = 0; ij < L; ij++ ){
            size_t y = kmer[ij] % Global::A2powerK[K_+1];
            for (size_t j = 0; j < W_; j++) {
                // parallelize for: n_[K_][y][j] += r_[n][L-ij+j-1];
                atomic_float_add( &(n_[K_][y][j]), r_[n][L+padding_-ij+j-1] );
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

    // update model parameters v[k][y][j] with updated k-mer counts,
    // alphas and model order
    motif_->updateV( n_, A_ );
}

void EM::updatePrior(){

    for( size_t n = 0; n < seqs_.size(); n++ ){
        size_t 	L = seqs_[n]->getL();
        size_t 	LW1 = L - W_ + 1;
        prior_[n][0] = 1.0f - q_;                         // probability of having no motif on the sequence
        float pos_i = q_ / static_cast<float>( LW1 );   // probability of having a motif at position i on the sequence
        for( size_t i = 1; i <= L; i++ ){
            prior_[n][i] = pos_i;
        }
    }

}

void EM::EStep_slow(){
    /**
     * together with MStep_slow(), this is the slow version of EM
     */

    float llikelihood = 0.f;

    motif_->calculateLinearS( bgModel_->getV() );

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
            r_[n][i] *= prior_[n][i];
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

}

void EM::MStep_slow(){
    /**
     * together with EStep_slow(), this is the slow version of EM
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

int EM::mask(){
    /**
     * upgraded version of EM to eliminate the effect of unrelated motifs
     */
    auto t0_wall = std::chrono::high_resolution_clock::now();

    /**
     * E-step for k = 0 to estimate weights r
     */
    // calculate the log odds ratio for k=0
    for( size_t y = 0; y < Global::A2powerK[1]; y++ ){
        for( size_t j = 0; j < W_; j++ ){
            s_[y][j] = motif_->getV()[0][y][j]
                       / bgModel_->getV()[0][y];
        }
    }

    // initialized positional priors
    updatePrior();
    std::vector<float> r_all;   // add all r's to an array for sorting
    size_t N_position = 0;      // count the number of all possible positions

    // calculate responsibilities r at all LW1 positions on sequence n
    // n runs over all sequences
//#pragma omp parallel for
    for( size_t n = 0; n < seqs_.size(); n++ ) {

        size_t  L = seqs_[n]->getL();
        size_t  LW1 = L - W_ + 1;
        size_t* kmer = seqs_[n]->getKmer();

        // initialize r_[n][i]
        for( size_t i = 1; i <= LW1; i++ ){
            r_[n][i] = 1.0f;
        }

        // when p(z_n > 0), it runs over all positions in sequence
        // in the region of [1, LW1]
        for( size_t i = 1; i <= LW1; i++ ){
            for( size_t j = 0; j < W_; j++ ){
                // extract (K+1)-mer y from positions (ij-K,...,ij)
                size_t y = kmer[i-1+j] % Global::A2powerK[1];
                r_[n][i] *= s_[y][j];
            }
        }

        // calculate the responsibilities and sum them up
        for( size_t i = 0; i <= LW1; i++ ){
            r_[n][i] *= prior_[n][i];
        }

        // add all un-normalized r to a vector
        for( size_t i = 0; i <= LW1; i++ ) {
            r_all.push_back( r_[n][i] );
        }

        // count the number of all possible motif positions
        N_position += LW1;
    }

    /**
     * found the cutoff of responsibility that covers the top f_% of
     * the motif occurrences
     * by applying partition-based selection algorithm, which runs in O(n)
     */

    // Sort log odds scores in descending order
    std::sort( r_all.begin(), r_all.end(), std::greater<float>() );

    // find the cutoff with f_% best r's
    float r_cutoff = r_all[size_t( N_position * Global::fraction )];
    // avoid r_cutoff from being too small
    r_cutoff = std::max(r_cutoff, Global::rCutoff);

    std::cout << "r cutoff: " << r_cutoff << std::endl;

    // create a vector to store all the r's which are above the cutoff
    std::vector<std::vector<size_t>> ridx( seqs_.size() );

    // put the index of the f_% r's in an array
    for( size_t n = 0; n < seqs_.size(); n++ ){
        for( size_t i = 1; i <= seqs_[n]->getL()-W_+1; i++ ){
            if( r_[n][i] >= r_cutoff ){
                ridx[n].push_back( i );
            } else {
                r_[n][i] = 0.f;
            }
        }
    }

    /**
     * optimize the model with the highest order from the top f_% global occurrences
     * using EM till convergence
     */
    bool 	iterate = true;
    size_t  iteration = 0;
    float 	v_diff, llikelihood_prev, llikelihood_diff;

    std::vector<std::vector<float>> v_before;
    v_before.resize( Global::A2powerK[K_+1] );
    for( size_t y = 0; y < Global::A2powerK[K_+1]; y++ ){
        v_before[y].resize( W_ );
    }

    // iterate over
    while( iterate && ( iteration < Global::maxIterations ) ){

        // get parameter variables with highest order before EM
        llikelihood_prev = llikelihood_;
        for( size_t y = 0; y < Global::A2powerK[K_+1]; y++ ){
            for( size_t j = 0; j < W_; j++ ){
                v_before[y][j] = motif_->getV()[K_][y][j];
            }
        }

        /**
         * E-step for f_% motif occurrences
         */
        float llikelihood = 0.0f;

        motif_->calculateLinearS( bgModel_->getV() );

        // calculate responsibilities r at all LW1 positions on sequence n
        // n runs over all sequences
#pragma omp parallel for reduction(+:llikelihood)
        for( size_t n = 0; n < seqs_.size(); n++ ){
            if(ridx[n].size() > 0) {
                size_t *kmer = seqs_[n]->getKmer();
                float normFactor = 0.f;

                // initialize all filtered r_[n][i]
                for (size_t idx = 0; idx < ridx[n].size(); idx++) {
                    r_[n][ridx[n][idx]] = 1.0f;
                }

                // update r based on the log odd ratios
                for (size_t idx = 0; idx < ridx[n].size(); idx++) {
                    for (size_t j = 0; j < W_; j++) {
                        size_t y = kmer[ridx[n][idx] + j - 1] % Global::A2powerK[K_ + 1];
                        r_[n][ridx[n][idx]] *= s_[y][j];
                    }
                    // calculate the responsibilities and sum them up
                    r_[n][ridx[n][idx]] *= prior_[n][ridx[n][idx]];
                    normFactor += r_[n][ridx[n][idx]];
                }

                // normalize responsibilities
                normFactor += r_[n][0];
                r_[n][0] /= normFactor;
                for (size_t idx = 0; idx < ridx[n].size(); idx++) {
                    r_[n][ridx[n][idx]] /= normFactor;
                }

                // calculate log likelihood over all sequences
                llikelihood += logf(normFactor);
            }
        }

        llikelihood_ = llikelihood;

        /**
         * M-step for f_% motif occurrences
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
        // using the f_% r's
        // n runs over all sequences
#pragma omp parallel for
        for( size_t n = 0; n < seqs_.size(); n++ ){
            if(ridx[n].size() > 0) {
                size_t *kmer = seqs_[n]->getKmer();
                for (size_t idx = 0; idx < ridx[n].size(); idx++) {
                    for (size_t j = 0; j < W_; j++) {
                        size_t y = kmer[ridx[n][idx] + j - 1] % Global::A2powerK[K_ + 1];
                        // n_[K_][y][j] += r_[n][ridx[n][idx]];
                        atomic_float_add(&(n_[K_][y][j]), r_[n][ridx[n][idx]]);
                    }
                }
            }
        }

        // compute fractional occurrence counts from higher to lower order
        // k runs over all orders
        for( size_t k = K_; k > 0; k-- ){
            for( size_t y = 0; y < Global::A2powerK[k+1]; y++ ){
                size_t y2 = y % Global::A2powerK[k];
                for( size_t j = 0; j < W_; j++ ){
                    n_[k-1][y2][j] += n_[k][y][j];
                }
            }
        }

        // update model parameters v[k][y][j], due to the k-mer counts, alphas and model order
        motif_->updateV( n_, A_ );

        /**
         * check parameter difference for convergence
         */
        v_diff = 0.0f;
        for( size_t y = 0; y < Global::A2powerK[K_+1]; y++ ){
            for( size_t j = 0; j < W_; j++ ){
                v_diff += fabsf( motif_->getV()[K_][y][j] - v_before[y][j] );
            }
        }

        // check the change of likelihood for convergence
        llikelihood_diff = llikelihood_ - llikelihood_prev;
        if( Global::verbose ) {
            std::cout << "i=" << iteration
                      << std::fixed
                      << "\tllh=" << llikelihood_
                      << "\tllh_diff=" << llikelihood_diff
                      << "\t\tmodel_diff=" << v_diff
                      << std::endl;
        }
        if( v_diff < Global::EMepsilon )				iterate = false;
        if( llikelihood_diff < 0 and iteration > 10 )	iterate = false;

        iteration++;
    }

    // calculate probabilities
    motif_->calculateP();

    auto t1_wall = std::chrono::high_resolution_clock::now();
    auto t_diff = std::chrono::duration_cast<std::chrono::duration<double>>(t1_wall-t0_wall);
    std::cout << "\n--- Runtime for EM: " << t_diff.count()
              << " seconds ("<< std::to_string(iteration)
              << " iterations) ---\n";

    return 0;
}

void EM::optimizeQ(){

    float N_0 = 0.f;
    size_t posN = seqs_.size();
    for (size_t n = 0; n < posN; n++) {
        N_0 += r_[n][0];
    }

    std::cout << "N_0 = " << N_0 << std::endl;

    q_ = ( posN - N_0 + 1.f ) / ( posN + 2.f );

}


void EM::optimizePos() {

    // todo: note this function currently only works for sequences of the same length

    size_t posN = seqs_.size();
    size_t LW1_ = seqs_[0]->getL()-W_+1;
    // update positional prior using the smoothness parameter beta
    // according to Eq. 149
    // calculate the value N0 and Ni_i:
    N0_ = 0.f;
    for (size_t n = 0; n < posN; n++) {
        N0_ += r_[n][0];
    }

    for( size_t i = 0; i < LW1_; i++ ) {
        // Ni_[i] = 0;
        for (size_t n = 0; n < posN; n++) {
            Ni_[i] += r_[n][seqs_[n]->getL()-i];
        }
    }

    size_t method_flag = 4;

    if( method_flag == 1 ){

        // method 1: Using a flat Bayesian prior on positional preference
        // according to Eq. 149
        for( size_t i = 0; i < LW1_; i++ ) {
            // update pi according to Eq. 149
            pi_[i+1] = (Ni_[i] + beta1_ - 1.f) / (posN - N0_ + LW1_ * (beta1_ - 1.f));
        }

    } else if ( method_flag == 2 ){

        // method 2: Using a prior for penalising jumps in the positional
        // preference profile
        // run a few iterations of conjugate gradients (e.g. 5~10)
        for( size_t i = 0; i < LW1_; i++ ) {
            b_vector_[i] = ( Ni_[i] - ( posN - N0_ ) * pi_[i+1] ) / beta2_;
        }

        Eigen::ConjugateGradient<Eigen::MatrixXf, Eigen::Lower | Eigen::Upper> cg;

        cg.compute( A_matrix_ );
        cg.setMaxIterations( 5 );
        si_ = cg.solve( b_vector_ );

        // normalize pi
        norm_ = 1.e-5f;
        for( size_t i = 1; i <= LW1_; i++ ) {
            pi_[i] = expf( si_[i-1] );
            norm_ += pi_[i];
        }
        for( size_t i = 1; i <= LW1_; i++ ) {
            pi_[i] /= norm_;
        }

        // calculate objective function:
        float Q;
        float p1 = 0.f;
        float p2 = 0.f;
        float sumexp = 0.f;
        // calculate objective function
        for (size_t i = 0; i < LW1_; i++) {
            sumexp += expf(si_[i]);
        }
        for (size_t i = 0; i < LW1_; i++) {
            p1 += Ni_[i] * (si_[i] - logf(sumexp));
        }
        for (size_t i = 2; i < LW1_; i++) {
            p2 += beta2_ * powf(si_[i] - si_[i - 1], 2.f);
        }

        Q = p1 - p2;

        if( Global::verbose ){
            std::cout << "cost function = " << Q << std::endl;
        }
        // ======== end of objfun calculation =======

        // method 2b: update smoothness parameter beta using positional prior distribution
        // from the data
        // according to Eq. 158, note the indices for si
        float sum = 1.e-8f;
        for (size_t i = 1; i < LW1_; i++) {
            sum += powf( si_[i] - si_[i-1], 2.f );
        }

        // a) update beta by its expectation value
        //beta2_ = LW1_+1 / sum;

        // b) sample beta2 from the Gamma distribution
        std::default_random_engine generator;
        std::gamma_distribution<double> distribution( (LW1_+1)/2, sum/2);
        beta2_ = distribution(generator);

        if( Global::verbose ) {
            std::cout << "N0=" << N0_ << ", N=" << posN
                      << ", LW1=" << LW1_ << ", norm=" << norm_
                      << ", sum=" << sum << ", beta=" << (LW1_ + 1) / sum
                      << std::endl;
        }

    } else if( method_flag == 3 ) {
        // prior penalising kinks in the positional preference profile
        for (size_t i = 0; i < LW1_; i++) {
            b_vector_[i] = (Ni_[i] - (posN - N0_) * pi_[i+1]) * 4 / beta2_;
        }

        Eigen::ConjugateGradient<Eigen::MatrixXf, Eigen::Lower | Eigen::Upper> cg;
        Eigen::MatrixXf B_matrix = getBmatrix(LW1_);
        cg.compute(B_matrix);
        cg.setMaxIterations(5);
        si_ = cg.solve(b_vector_);

        // normalize pi
        norm_ = 1.e-5f;
        for (size_t i = 1; i <= LW1_; i++) {
            pi_[i] = expf(si_[i-1]);
            norm_ += pi_[i];
        }
        for (size_t i = 1; i <= LW1_; i++) {
            pi_[i] /= norm_;
        }

        float sum = 1.e-5f;
        for (size_t i = 1; i < LW1_-1; i++) {
            sum += powf(si_[i] - (si_[i-1] + si_[i+1]) / 2, 2.f);
        }
        // a) update beta by its expectation value 2
        beta2_ = LW1_ / sum;

        // b) update beta by its expectation value 2
        //beta2_ = ( LW1-2 ) / sum;

        if( Global::verbose ) {
            std::cout << "N0=" << N0_ << ", N=" << posN
                      << ", LW1=" << LW1_ << ", norm=" << norm_
                      << ", sum=" << sum << ", beta=" << (LW1_ + 1) / sum
                      << std::endl;
        }
    } else if( method_flag == 4 ){

        // Use L-BFGS as optimizer

        // pre-define parameters
        LBFGSpp::LBFGSParam<float> param;
        param.epsilon   = 1.e-2f;
        param.delta     = 1.e-4f;
        param.past      = 1;
        param.max_iterations = 1000;

        // create colver and function object
        LBFGSpp::LBFGSSolver<float> solver( param );

        // solve the function
        Eigen::VectorXf grad = Eigen::VectorXf::Zero( LW1_ );

        // set up the objective function
        ObjFun func( LW1_, posN, N0_, beta2_, Ni_, A_matrix_ );

        float fx;
        int niter;
        niter = solver.minimize( func, si_, fx );
        if( Global::verbose ){

        }

        // calculate beta2 after optimization of si
        float sumS = 0.f;
        for( size_t i = 1; i < LW1_; i++){
            sumS += powf(si_[i] - si_[i-1], 2.f);
        }
        beta2_ = ( LW1_-1 ) / sumS;

        if( Global::verbose ){
            std::cout << niter << " iterations, f(x)=" << fx
                      << ", optimized beta2 = " << beta2_ << std::endl;
        }
        // update pi
        norm_ = 1.e-5f;
        for (size_t i = 1; i <= LW1_; i++) {
            pi_[i] = expf(si_[i-1]);
            norm_ += pi_[i];
        }
        for (size_t i = 1; i <= LW1_; i++) {
            pi_[i] /= norm_;
        }

    }

    // update pos_ni by pi_i
    for( size_t i = 1; i <= LW1_; i++ ) {
        for (size_t n = 0; n < posN; n++) {
            prior_[n][i] = pi_[i] * q_;
        }
    }

}


float EM::obj_fun( Eigen::VectorXf& si, Eigen::VectorXf& grad ){

    float Q;
    float p1 = 0.f;
    float p2 = 0.f;
    float sumexp = 0.f;
    size_t LW1_ = seqs_[0]->getL()-W_+1;

    // calculate objective function
    for( size_t i = 0; i < LW1_; i++ ){
        sumexp += expf( si[i] );
    }
    for( size_t i = 0; i < LW1_; i++ ){
        p1 += Ni_[i] * ( si[i] - logf( sumexp ) );
    }
    for( size_t i = 2; i < LW1_; i++ ){
        p2 += beta2_ * powf( si[i] - si[i-1], 2.f );
    }

    Q = p2 - p1;

    // update the gradients
    Eigen::VectorXf dot_product = A_matrix_ * si;
    for( size_t i = 0; i < LW1_; i++ ){
        grad[i] = - ( Ni_[i] - (seqs_.size() - N0_) * expf( si[i] )
                               / sumexp - beta2_ * dot_product[i] );
    }

    return Q;
}

Eigen::MatrixXf EM::getAmatrix( size_t w ) {
    // define matrix A:
    Eigen::MatrixXf A_matrix = Eigen::MatrixXf::Zero( w, w );
    for( size_t i = 0; i < w; i++ ){
        for( size_t j = 0; j < w; j++ ){
            if( i == j ) {
                A_matrix(i, j) = 2;
            } else if( abs( i-j ) == 1 ){
                A_matrix(i, j) = -1;
            }
        }
    }
    A_matrix(0, 0) = 1;
    A_matrix(w-1, w-1) = 1;

    return A_matrix;
}

Eigen::MatrixXf EM::getBmatrix( size_t w ) {
    // define matrix B:
    Eigen::MatrixXf B_matrix = Eigen::MatrixXf::Zero( w, w );

    for( size_t i = 0; i < w; i++ ){
        for( size_t j = 0; j < w; j++ ){
            if( i == j ) {
                B_matrix(i, j) = 6;
            } else if( std::abs( i-j ) == 1 ){
                B_matrix(i, j) = -4;
            } else if( std::abs( i-j ) == 2 ){
                B_matrix(i, j) = 1;
            }
        }
    }
    B_matrix(0, 0) = B_matrix(w-1, w-1) = 1;
    B_matrix(1, 1) = B_matrix(w-2, w-2) = 5;
    B_matrix(0, 1) = B_matrix(1, 0) = B_matrix(w-1, w-2) = B_matrix(w-2, w-1) = -2;

    return B_matrix;
}

float* EM::getPi(){
    return pi_;
}


float** EM::getR(){
    return r_;
}

float EM::getQ(){
    return q_;
}



void EM::print(){

    std::cout << std::fixed << std::setprecision( 2 );

    std::cout << " _________________________" << std::endl
              << "|                         |" << std::endl
              << "| k-mer Fractional Counts |" << std::endl
              << "|_________________________|" << std::endl
              << std::endl;

    for( size_t j = 0; j < W_; j++ ){
        for( size_t k = 0; k < K_+1; k++ ){
            float sum = 0.f;
            for( size_t y = 0; y < Global::A2powerK[k+1]; y++ ) {
                sum += n_[k][y][j];
                std::cout << n_[k][y][j] << '\t';
            }
            std::cout << "sum=" << sum << std::endl;
        }
    }

    std::cout << std::fixed << std::setprecision( 4 );

    std::cout << " ___________________________" << std::endl
              << "|                           |" << std::endl
              << "| Conditional Probabilities |" << std::endl
              << "|___________________________|" << std::endl
              << std::endl;

    for( size_t j = 0; j < W_; j++ ){
        for( size_t k = 0; k < K_+1; k++ ){
            float sum = 0.f;
            for( size_t y = 0; y < Global::A2powerK[k+1]; y++ ) {
                sum += motif_->getV()[k][y][j];
                std::cout << motif_->getV()[k][y][j] << '\t';
            }
            std::cout << "sum=" << sum << std::endl;
//            assert( fabsf( sum - Global::A2powerK[k] ) < 1.e-3f );
        }
    }

    std::cout << " ____________________" << std::endl
              << "|                    |" << std::endl
              << "| Full Probabilities |" << std::endl
              << "|____________________|" << std::endl
              << std::endl;

    for( size_t j = 0; j < W_; j++ ){
        for( size_t k = 0; k < K_+1; k++ ){
            float sum = 0.f;
            for( size_t y = 0; y < Global::A2powerK[k+1]; y++ ) {
                sum += motif_->getP()[k][y][j];
                std::cout << motif_->getP()[k][y][j] << '\t';
            }
            std::cout << "sum=" << sum << std::endl;
//            assert( fabsf( sum-1.0f ) < 1.e-4f );
        }
    }

    std::cout << " _______________" << std::endl
              << "|               |" << std::endl
              << "| Motifs Scores |" << std::endl
              << "|_______________|" << std::endl
              << std::endl;

    for( size_t j = 0; j < W_; j++ ){
        float sum = 0.f;
        for( size_t y = 0; y < Global::A2powerK[K_+1]; y++ ) {
            sum += s_[y][j];
            std::cout << s_[y][j] << '\t';
        }
        std::cout << "sum=" << sum << std::endl;
    }

    std::cout << std::endl << "Fraction parameter q="
              << std::setprecision( 4 ) << q_ << std::endl;
}

void EM::printR(){

    // print out weights r for each position on each sequence
    std::cout << " _____________" << std::endl
              << "|             |" << std::endl
              << "|  Posterior  |" << std::endl
              << "|_____________|" << std::endl
              << std::endl;

    if( !Global::slowEM ){
        for( size_t n = 0; n < seqs_.size(); n++ ){
            std::cout << "seq " << n+1 << ":" << std::endl;
            size_t L = seqs_[n]->getL();
            float sum = 0.f;
            std::cout << r_[n][0] << '\t';
            sum += r_[n][0];
            for( size_t i = 1; i <= L-W_+1; i++ ){
                sum += r_[n][L+padding_-i];
                std::cout << r_[n][L+padding_-i] << '\t';
            }
            std::cout << /*"sum=" << sum << */std::endl;
            assert( fabsf( sum-1.0f ) < 1.e-4f );
        }
    } else {
        for( size_t n = 0; n < seqs_.size(); n++ ){
            std::cout << "seq " << n+1 << ":" << std::endl;
            size_t L = seqs_[n]->getL();
            float sum = 0.f;
            for( size_t i = 0; i <= L-W_+1; i++ ){
                sum += r_[n][i];
                std::cout << r_[n][i] << '\t';
            }
            std::cout << /*"sum=" << sum <<*/ std::endl;
            assert( fabsf( sum-1.0f ) < 1.e-4f );
        }
    }
}

void EM::write( char* odir, std::string basename ){

	/**
	 * 	 * save BaMM (hyper-)parameters in flat files:
	 * (1) posSequenceBasename.counts:	    refined fractional counts of (k+1)-mers
	 * (2) posSequenceBasename.weights:     responsibilities, posterior distributions
	 * (3) posSequenceBasename.positions:   position of motif(s) on each sequence
	 */

	std::string opath = std::string( odir ) + '/' + basename;

	// output (k+1)-mer counts n[k][y][j]
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

	// output position(s) of motif(s): prior_[n][i]
	std::string opath_pos = opath + ".positions";
	std::ofstream ofile_pos( opath_pos.c_str() );

	ofile_pos << "seq\tlength\tstrand\tstart..end\tpattern" << std::endl;

    float cutoff = 0.3f;	// threshold for having a motif on the sequence
                            // in terms of responsibilities
    for( size_t n = 0; n < seqs_.size(); n++ ){

        size_t L = seqs_[n]->getL();
        size_t raw_L = Global::ss ? L : ( L-1 ) / 2;

        for( size_t i = 0; i < L-W_+1; i++ ){
            if( r_[n][L+padding_+1-i] >= cutoff ){
                ofile_pos << seqs_[n]->getHeader() << '\t' << raw_L << '\t'
                          << ( ( i < raw_L ) ? '+' : '-' ) << '\t' << i+1 << ".." << i+W_ << '\t';
                for( size_t b = i; b < i+W_; b++ ){
                    ofile_pos << Alphabet::getBase( seqs_[n]->getSequence()[b] );
                }
                ofile_pos << std::endl;
            }
        }
    }
}