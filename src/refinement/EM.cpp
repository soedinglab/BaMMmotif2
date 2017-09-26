//
// Created by wanwan on 16.08.17.
//
#include "EM.h"

EM::EM( Motif* motif, BackgroundModel* bgModel, std::vector<Sequence*> seqs, float q, bool optimizeQ ){

    motif_      = motif;
	bgModel_    = bgModel;
	q_          = q;
	seqs_       = seqs;
    optimizeQ_ = optimizeQ;
	// get motif (hyper-)parameters from motif class
	K_ = motif_->getK();
	W_ = motif_->getW();
	Y_ = motif_->getY();
	s_ = motif_->getS();
	A_ = motif_->getA();
	K_bg_ = ( bgModel_->getOrder() < K_ ) ? bgModel_->getOrder() : K_;

	// allocate memory for r_[n][i], pos_[n][i], z_[n]
	r_ = ( float** )calloc( seqs_.size(), sizeof( float* ) );
	pos_ = ( float** )calloc( seqs_.size(), sizeof( float* ) );
	for( size_t n = 0; n < seqs_.size(); n++ ){
		size_t LW2 = seqs_[n]->getL() - W_ + 2;
		r_[n] = ( float* )calloc( LW2, sizeof( float ) );
		pos_[n] = ( float* )calloc( LW2, sizeof( float ) );
	}

	// allocate memory for n_[k][y][j]
	n_ = ( float*** )calloc( K_+1, sizeof( float** ) );
	for( size_t k = 0; k < K_+1; k++ ){
		n_[k] = ( float** )calloc( Y_[k+1], sizeof( float* ) );
		for( size_t y = 0; y < Y_[k+1]; y++ ){
			n_[k][y] = ( float* )calloc( W_, sizeof( float ) );
		}
	}
}

EM::~EM(){
    for( size_t n = 0; n < seqs_.size(); n++ ){
        free( r_[n] );
        free( pos_[n] );
    }
    free( r_ );
    free( pos_ );

    for( size_t k = 0; k < K_+1; k++ ){
        for( size_t y = 0; y < Y_[k+1]; y++ ){
            free( n_[k][y] );
        }
        free( n_[k] );
    }
    free( n_ );
}

int EM::optimize(){

    clock_t t0 = clock();

    bool 	iterate = true;

    float 	v_diff, llikelihood_prev, llikelihood_diff;

    std::vector<std::vector<float>> v_before;
    v_before.resize( Y_[K_+1] );
    for( size_t y = 0; y < Y_[K_+1]; y++ ){
        v_before[y].resize( W_ );
    }

    size_t iteration = 0;

    // iterate over
    while( iterate && ( iteration < maxEMIterations_ ) ){

        iteration++;

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

        // optimize hyperparameter q at every 10 steps
        if( optimizeQ_ and iteration % 10 == 1 )    optimize_q();

        // check parameter difference for convergence
        v_diff = 0.0f;
        for( size_t y = 0; y < Y_[K_+1]; y++ ){
            for( size_t j = 0; j < W_; j++ ){
                v_diff += fabsf( motif_->getV()[K_][y][j] - v_before[y][j] );
            }
        }

        // check the change of likelihood for convergence
        llikelihood_diff = llikelihood_ - llikelihood_prev;

        std::cout << iteration << "th iteration, delta_llikelihood=" << llikelihood_diff << std::endl;

        if( v_diff < epsilon_ )							iterate = false;
        if( !optimizeQ_ and llikelihood_diff < 0 and iteration > 10 )	iterate = false;

/*        // for making a movie out of all iterations
        // calculate probabilities
        motif_->calculateP();
        // write out the learned model
        motif_->write( odir, "ChIP_GATA3_iter_" + std::to_string( iteration ) );*/
    }

    // calculate probabilities
    motif_->calculateP();

    fprintf( stdout, "\n--- Runtime for EM: %.4f seconds ---\n",
             ( ( float )( clock() - t0 ) ) / CLOCKS_PER_SEC );
    return 0;
}

void EM::EStep(){

    llikelihood_ = 0.0f;

    motif_->calculateLinearS( bgModel_->getV(), K_bg_ );

    // todo: parallelize the code
//	#pragma omp parallel for

    // calculate responsibilities r at all LW1 positions on sequence n
    // n runs over all sequences
    for( size_t n = 0; n < seqs_.size(); n++ ){

        size_t 	L = seqs_[n]->getL();
        size_t 	LW1 = L - W_ + 1;
        size_t*	kmer = seqs_[n]->getKmer();
        float 	normFactor = 0.0f;

        // initialize r_[n][i] and pos_[n][i]
        float pos_i = q_ / static_cast<float>( LW1 );
        for( size_t i = 1; i <= LW1; i++ ){
            r_[n][i] = 1.0f;
            pos_[n][i] = pos_i;
        }
        pos_[n][0] = 1.0f - q_;
        r_[n][0] = pos_[n][0];

        // when p(z_n > 0), ij = i+j runs over all positions in sequence
        for( size_t ij = 0; ij < L; ij++ ){

            // extract (K+1)-mer y from positions (i-k,...,i)
            size_t y = kmer[ij] % Y_[K_+1];
            // j runs over all motif positions
            size_t padding = ( static_cast<int>( ij-L+W_ ) > 0 ) * ( ij-L+W_ );
            for( size_t j = padding; j < ( W_ < (ij+1) ? W_ : ij+1 ); j++ ){
                r_[n][LW1-ij+j] *= s_[y][j];
            }
        }

        // calculate the responsibilities and sum them up
        normFactor += r_[n][0];
        for( size_t i = 1; i <= LW1; i++ ){
            r_[n][i] *= pos_[n][LW1+1-i];
            normFactor += r_[n][i];
        }

        // normalize responsibilities
        for( size_t i = 0; i <= LW1; i++ ){
            r_[n][i] /= normFactor;
        }

        // calculate log likelihood over all sequences
        llikelihood_ += logf( normFactor );
    }

}

void EM::MStep(){

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
    for( size_t n = 0; n < seqs_.size(); n++ ){
        size_t L = seqs_[n]->getL();
        size_t* kmer = seqs_[n]->getKmer();

        // ij = i+j runs over all positions i on sequence n
        for( size_t ij = 0; ij < L; ij++ ){

            size_t y = kmer[ij] % Y_[K_+1];

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

    // update model parameters v[k][y][j], due to the kmer counts, alphas and model order
    motif_->updateV( n_, A_, K_ );
}

int EM::advance() {
    /**
     * upgraded version of EM to eliminate the effect of unrelated motifs
     */
    clock_t t0 = clock();

    /**
     * run 5 steps of EM for k=0 and optimize q in the meantime
     */
    for( size_t i = 0; i < 1; i++ ){

        llikelihood_ = 0.0f;

        /**
         * E-step for k = 0
         */
        // calculate the log odds ratio for k=0
        for( size_t y = 0; y < Y_[1]; y++ ){
            for( size_t j = 0; j < W_; j++ ){
                s_[y][j] = motif_->getV()[0][y][j] / bgModel_->getV()[0][y];
            }
        }

        // todo: parallelize the code
        // calculate responsibilities r at all LW1 positions on sequence n
        // n runs over all sequences
        for( size_t n = 0; n < seqs_.size(); n++ ) {

            size_t L = seqs_[n]->getL();
            size_t LW1 = L - W_ + 1;
            size_t *kmer = seqs_[n]->getKmer();
            float normFactor = 0.0f;

            // initialize r_[n][i] and pos_[n][i]
            float pos_i = q_ / static_cast<float>( LW1 );
            for (size_t i = 1; i <= LW1; i++) {
                r_[n][i] = 1.0f;
                pos_[n][i] = pos_i;
            }
            pos_[n][0] = 1.0f - q_;
            r_[n][0] = pos_[n][0];

            // when p(z_n > 0), ij = i+j runs over all positions in sequence
            for (size_t ij = 0; ij < L; ij++) {

                // extract monomer y at position i
                size_t y = kmer[ij] % Y_[1];
                // j runs over all motif positions
                size_t padding = (static_cast<int>( ij - L + W_ ) > 0) * (ij - L + W_);
                for (size_t j = padding; j < (W_ < (ij + 1) ? W_ : ij + 1); j++) {
                    r_[n][LW1 - ij + j] *= s_[y][j];
                }
            }

            // calculate the responsibilities and sum them up
            normFactor += r_[n][0];
            for (size_t i = 1; i <= LW1; i++) {
                r_[n][i] *= pos_[n][LW1 + 1 - i];
                normFactor += r_[n][i];
            }

            // normalize responsibilities
            for (size_t i = 0; i <= LW1; i++) {
                r_[n][i] /= normFactor;
            }

            // calculate log likelihood over all sequences
            llikelihood_ += logf(normFactor);
        }

        /**
         * M-Step for k=0
         */
        // reset the fractional counts n for k = 0
        for( size_t y = 0; y < Y_[1]; y++ ){
            for( size_t j = 0; j < W_; j++ ){
                n_[0][y][j] = 0.0f;
            }
        }

        // compute fractional occurrence counts for k = 0
        // n runs over all sequences
        for( size_t n = 0; n < seqs_.size(); n++ ){
            size_t L = seqs_[n]->getL();
            size_t* kmer = seqs_[n]->getKmer();

            // ij = i+j runs over all positions i on sequence n
            for( size_t ij = 0; ij < L; ij++ ){
                size_t y = kmer[ij] % Y_[1];
                size_t padding = ( static_cast<int>( ij-L+W_ ) > 0 ) * ( ij-L+W_ );
                for( size_t j = padding; j < ( W_ < (ij+1) ? W_ : ij+1 ); j++ ){
                    n_[0][y][j] += r_[n][L-W_+1-ij+j];
                }
            }
        }

        // update model parameters v[k][y][j], due to the kmer counts, alphas and model order
        motif_->updateV( n_, A_, 0 );

        /**
         * optimize fractional prior q
         */
        optimize_q();
        // check the optimized q
        std::cout << "optimized q = " << q_ << std::endl;
    }

    /**
     * found the cutoff of responsibility that covers the top 10% of the motif occurrences
     * by applying partion-based selection algorithm, which runs in O(n)
     */
    std::vector<float> r_all;
    size_t pos_count = 0;
    // add all r's to an array for sorting
    for( size_t n = 0; n < seqs_.size(); n++ ){
        size_t LW1 = seqs_[n]->getL() - W_ + 1;
        for (size_t i = 1; i <= LW1; i++) {
            r_all.push_back( r_[n][i] );
        }
        pos_count += LW1;
    }
    // Sort log odds scores in descending order
    std::sort( r_all.begin(), r_all.end(), std::greater<float>() );
    float r_cutoff = r_all[pos_count / 20];
    std::vector<std::vector<size_t>> ri;
    ri.resize( seqs_.size() );

    for( size_t n = 0; n < seqs_.size(); n++ ){
        size_t LW1 = seqs_[n]->getL() - W_ + 1;
        r_all.push_back( r_[n][0] );
        for( size_t i = 1; i <= LW1; i++ ){
            if( r_[n][i] >= r_cutoff ){
                ri[n].push_back( i );
            }
        }
    }

    /**
     * optimize the model from the top 10% global occurrences with higher orders till convergence
     */
    bool 	iterate = true;

    float 	v_diff, llikelihood_prev, llikelihood_diff;

    std::vector<std::vector<float>> v_before;
    v_before.resize( Y_[K_+1] );
    for( size_t y = 0; y < Y_[K_+1]; y++ ){
        v_before[y].resize( W_ );
    }

    size_t iteration = 0;

    // iterate over
    while( iterate && ( iteration < maxEMIterations_ ) ){

        iteration++;

        // get parameter variables with highest order before EM
        llikelihood_prev = llikelihood_;
        for( size_t y = 0; y < Y_[K_+1]; y++ ){
            for( size_t j = 0; j < W_; j++ ){
                v_before[y][j] = motif_->getV()[K_][y][j];
            }
        }

        /**
         * E-step for 10% motif occurrences
         */
        llikelihood_ = 0.0f;

        motif_->calculateLinearS( bgModel_->getV(), K_bg_ );

        // todo: parallelize the code
    //	#pragma omp parallel for

        // calculate responsibilities r at all LW1 positions on sequence n
        // n runs over all sequences
        for( size_t n = 0; n < seqs_.size(); n++ ){

            size_t 	L = seqs_[n]->getL();
            size_t 	LW1 = L - W_ + 1;
            size_t  LW2 = L - W_ + 2;
            size_t*	kmer = seqs_[n]->getKmer();
            float 	normFactor = 0.0f;

            // initialize r_[n][i] and pos_[n][i]
            float pos_i = q_ / static_cast<float>( LW1 );
            for( size_t idx = 0; idx < ri[n].size(); idx++ ){
                r_[n][ri[n][idx]] = 1.0f;
                pos_[n][ri[n][idx]] = pos_i;
            }
            pos_[n][0] = 1.0f - q_;
            r_[n][0] = pos_[n][0];


            for( size_t idx = 0; idx < ri[n].size(); idx++ ){
                for( size_t j = 0; j < W_; j++ ){
                    size_t y = kmer[LW2-ri[n][idx]+j] % Y_[K_+1];
                    r_[n][ri[n][idx]] *= s_[y][j];
                }
            }

            // calculate the responsibilities and sum them up
            normFactor += r_[n][0];

            for( size_t idx = 0; idx < ri[n].size(); idx++ ){
                r_[n][ri[n][idx]] *= pos_[n][LW2-ri[n][idx]];
                normFactor += r_[n][ri[n][idx]];
            }
            // normalize responsibilities
            r_[n][0] /= normFactor;
            for( size_t idx = 0; idx < ri[n].size(); idx++ ){
                r_[n][ri[n][idx]] /= normFactor;
            }

            // calculate log likelihood over all sequences
            llikelihood_ += logf( normFactor );
        }

        /**
         * M-step for 10% motif occurrences
         */
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
        for( size_t n = 0; n < seqs_.size(); n++ ){
            size_t  L = seqs_[n]->getL();
            size_t  LW2 = L - W_ + 2;
            size_t* kmer = seqs_[n]->getKmer();

            for( size_t idx = 0; idx < ri[n].size(); idx++ ){
                for( size_t j = 0; j < W_; j++ ){
                    size_t y = kmer[LW2-ri[n][idx]+j] % Y_[K_+1];
                    n_[K_][y][j] += r_[n][ri[n][idx]];
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

        // update model parameters v[k][y][j], due to the kmer counts, alphas and model order
        motif_->updateV( n_, A_, K_ );

        /**
         * check parameter difference for convergence
         */
        v_diff = 0.0f;
        for( size_t y = 0; y < Y_[K_+1]; y++ ){
            for( size_t j = 0; j < W_; j++ ){
                v_diff += fabsf( motif_->getV()[K_][y][j] - v_before[y][j] );
            }
        }

        // check the change of likelihood for convergence
        llikelihood_diff = llikelihood_ - llikelihood_prev;
        std::cout << iteration << "th iteration, delta_llikelihood=" << llikelihood_diff << std::endl;
        if( v_diff < epsilon_ )							iterate = false;
        if( llikelihood_diff < 0 and iteration > 10 )	iterate = false;

    }

    // calculate probabilities
    motif_->calculateP();

    fprintf( stdout, "\n--- Runtime for EM: %.4f seconds ---\n",
             ( ( float )( clock() - t0 ) ) / CLOCKS_PER_SEC );
    return 0;
}
void EM::optimize_q(){

    float N0 = 0;                   // expectation value of the count of sequences without a query motif

    for( size_t n = 0; n < seqs_.size(); n++ ){
        N0 += r_[n][0];
    }

    q_ = ( N0 + 1.f ) / ( ( float )seqs_.size() + 2.f );

}

float** EM::getR(){
    return r_;
}

float EM::getQ(){
    return q_;
}

void EM::print(){

    // print out motif parameter v
    for( size_t j = 0; j < W_; j++ ){
        for( size_t k = 0; k < K_+1; k++ ){
            for( size_t y = 0; y < Y_[k+1]; y++ ){
                std::cout << std::setprecision( 5 ) << motif_->getV()[k][y][j] << '\t';
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

}

void EM::write( char* odir, std::string basename, bool ss ){

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

	ofile_pos << "seq\tlength\tstrand\tstart..end\tpattern" << std::endl;

    float cutoff = 0.3f;	// threshold for having a motif on the sequence
                            // in terms of responsibilities
    for( size_t n = 0; n < seqs_.size(); n++ ){

        size_t L = seqs_[n]->getL();
        size_t LW1 = seqs_[n]->getL() - W_ + 1;

        if( ss ){

            for( size_t i = LW1; i > 0; i-- ){
                if( r_[n][i] >= cutoff ){
                    ofile_pos << '>' << seqs_[n]->getHeader() << '\t'
                              << L << '\t'
                              << '+' << '\t'
                              << LW1-i+1 << ".." << LW1-i+W_<< '\t';
                    for( size_t b = 0; b < W_; b++ ){
                        ofile_pos << Alphabet::getBase( seqs_[n]->getSequence()[LW1-i+b] );
                    }
                    ofile_pos << std::endl;
                }
            }

        } else {
            L = ( L - 1 ) / 2;
            for( size_t i = LW1; i > 0; i-- ){
                if( r_[n][i] >= cutoff ){
                    ofile_pos << '>' << seqs_[n]->getHeader() << '\t'
                              << L << '\t'
                              << ( ( i < L ) ? '-' : '+' ) << '\t'
                              << LW1-i+1 << ".." << LW1-i+W_<< '\t';
                    for( size_t b = 0; b < W_; b++ ){
                        ofile_pos << Alphabet::getBase( seqs_[n]->getSequence()[LW1-i+b] );
                    }
                    ofile_pos << std::endl;
                }
            }
        }

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