//
// Created by wanwan on 16.08.17.
//
#include "EM.h"

EM::EM( Motif* motif, BackgroundModel* bgModel,
        std::vector<Sequence*> seqs,
        float q, bool optimizeQ, bool verbose, float f ){

    motif_      = motif;
	bgModel_    = bgModel;
	q_          = q;
    f_          = f;
	seqs_       = seqs;
    optimizeQ_  = optimizeQ;

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
		r_[n] = ( float* )calloc( seqs_[n]->getL()/*-W_+1*/, sizeof( float ) );
		pos_[n] = ( float* )calloc( seqs_[n]->getL()/*-W_+1*/, sizeof( float ) );
	}

	// allocate memory for n_[k][y][j]
	n_ = ( float*** )calloc( K_+1, sizeof( float** ) );
	for( size_t k = 0; k < K_+1; k++ ){
		n_[k] = ( float** )calloc( Y_[k+1], sizeof( float* ) );
		for( size_t y = 0; y < Y_[k+1]; y++ ){
			n_[k][y] = ( float* )calloc( W_, sizeof( float ) );
		}
	}

    verbose_ = verbose;
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

        // optimize hyper-parameter q in the first 5 steps
        if( optimizeQ_ and iteration <= 5 )    optimize_q();

        // check parameter difference for convergence
        v_diff = 0.0f;
        for( size_t y = 0; y < Y_[K_+1]; y++ ){
            for( size_t j = 0; j < W_; j++ ){
                v_diff += fabsf( motif_->getV()[K_][y][j] - v_before[y][j] );
            }
        }

        // check the change of likelihood for convergence
        llikelihood_diff = llikelihood_ - llikelihood_prev;

        if( verbose_ ) std::cout << iteration << " iter, llh=" << llikelihood_
                                 << ", diff_llh=" << llikelihood_diff
                                 << ", v_diff=" << v_diff
                                 << std::endl;

        if( v_diff < epsilon_ )							iterate = false;
        if( llikelihood_diff < 0 and iteration > 10 )	iterate = false;

/*
        // for making a movie out of all iterations
        // calculate probabilities
        motif_->calculateP();
        // write out the learned model
        motif_->write( odir, filename + std::to_string( iteration ) );
*/

    }

    // calculate probabilities
    motif_->calculateP();

    fprintf( stdout, "\n--- Runtime for EM: %.4f seconds ---\n",
             ( ( float )( clock() - t0 ) ) / CLOCKS_PER_SEC );
    return 0;
}

void EM::EStep(){

    float llikelihood = 0.0f;

    motif_->calculateLinearS( bgModel_->getV(), K_bg_ );

    // calculate responsibilities r at all LW1 positions on sequence n
    // n runs over all sequences

#pragma omp parallel for reduction(+:llikelihood)
    for( size_t n = 0; n < seqs_.size(); n++ ){

        size_t 	L = seqs_[n]->getL();
        size_t 	LW1 = L - W_ + 1;
        size_t*	kmer = seqs_[n]->getKmer();
        float 	normFactor = 1.0f - q_;

        // initialize r_[n][i] and pos_[n][i]:
        // note here:   r_[n][0] and pos_[n][0] are for motif is absent
        //              and the indices of r_[n][i] is reverted, which means:
        //              r_[n][i] is the weight for motif being at position L-W-i on sequence n
        float pos_i = q_ / static_cast<float>( LW1 );
        for( size_t i = 0; i < LW1; i++ ){
            r_[n][i] = 1.0f;
            pos_[n][i] = pos_i;
        }

        // when p(z_n > 0), ij = i+j runs over all positions in sequence
        for( size_t ij = 0; ij < LW1; ij++ ){

            // extract (K+1)-mer y from positions (ij-K,...,ij)
            size_t y = kmer[ij] % Y_[K_+1];

            for( size_t j = 0; j < W_; j++ ){
                r_[n][L-W_-ij+j] *= s_[y][j];
            }

        }

        // calculate the responsibilities and sum them up
        for( size_t i = 0; i < LW1; i++ ){
            r_[n][i] *= pos_[n][i];
            normFactor += r_[n][i];
        }

        // normalize responsibilities
        for( size_t i = 0; i < LW1; i++ ){
            r_[n][i] /= normFactor;
        }

        // for the unused positions
        for( size_t i = LW1; i < L; i++ ){
            r_[n][i] = 0.0f;
       }
        // calculate log likelihood over all sequences
        llikelihood += logf( normFactor );
    }

    llikelihood_ = llikelihood;

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
//#pragma omp parallel for
    for( size_t n = 0; n < seqs_.size(); n++ ){
        size_t L = seqs_[n]->getL();
        size_t* kmer = seqs_[n]->getKmer();

        // ij = i+j runs over all positions i on sequence n
        for( size_t ij = 0; ij < L-W_+1; ij++ ){
            size_t y = kmer[ij] % Y_[K_+1];
            for( size_t j = 0; j < W_; j++ ){
//                __sync_fetch_and_add(n_[K_][y] + j, r_[n][L-W_-ij+j]);
                n_[K_][y][j] += r_[n][L-W_-ij+j];
            }
        }

    }

    // compute fractional occurrence counts from higher to lower order
    // k runs over all lower orders
    for( size_t k = K_; k > 0; k-- ){
        for( size_t y = 0; y < Y_[k+1]; y++ ){
            size_t y2 = y % Y_[k];
            for( size_t j = 0; j < W_; j++ ){
                n_[k-1][y2][j] += n_[k][y][j];
            }
        }
    }


    // update model parameters v[k][y][j] with updated kmer counts, alphas and model order
    motif_->updateV( n_, A_, K_ );
}

int EM::mask() {
    /**
     * upgraded version of EM to eliminate the effect of unrelated motifs
     */
    clock_t t0 = clock();

    /**
     * E-step for k = 0 to estimate weights r
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
        float normFactor = 1.0f - q_;

        // initialize r_[n][i] and pos_[n][i]
        float pos_i = q_ / static_cast<float>( LW1 );
        for (size_t i = 0; i < LW1; i++) {
            r_[n][i] = 1.0f;
            pos_[n][i] = pos_i;
        }


        // when p(z_n > 0), ij = i+j runs over all positions in sequence
        for (size_t ij = 0; ij < L; ij++) {

            // extract monomer y at position i
            size_t y = kmer[ij] % Y_[1];
            // j runs over all motif positions
            size_t padding = (static_cast<int>( ij - L + W_ ) > 0) * (ij - L + W_);
            for (size_t j = padding; j < (W_ < ij ? W_ : ij ); j++) {
                r_[n][L-W_-ij + j] *= s_[y][j];
            }
        }

        // calculate the responsibilities and sum them up
        for (size_t i = 0; i < LW1; i++) {
            r_[n][i] *= pos_[n][L-W_-i];
            normFactor += r_[n][i];
        }

        // normalize responsibilities
        for (size_t i = 0; i < LW1; i++) {
            r_[n][i] /= normFactor;
        }

        /**
         * optimize fractional prior q
         */
        if( optimizeQ_ )    optimize_q();

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
        for (size_t i = 0; i < LW1; i++) {
            r_all.push_back( r_[n][i] );
        }
        pos_count += LW1;
    }
    // Sort log odds scores in descending order
    std::sort( r_all.begin(), r_all.end(), std::greater<float>() );

    // find the cutoff with f_% best r's
    float r_cutoff = r_all[size_t((float)pos_count * f_)];

    std::vector<std::vector<size_t>> ri;
    ri.resize( seqs_.size() );

    // put the index of the f_% r's in an array
    for( size_t n = 0; n < seqs_.size(); n++ ){
        size_t LW1 = seqs_[n]->getL() - W_ + 1;
        for( size_t i = 0; i < LW1; i++ ){
            if( r_[n][i] >= r_cutoff ){
                ri[n].push_back( i );
            }
        }
    }

    /**
     * optimize the model from the top 10% global occurrences with the highest order till convergence
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
         * E-step for f_% motif occurrences
         */
        float llikelihood = 0.0f;

        motif_->calculateLinearS( bgModel_->getV(), K_bg_ );

        // calculate responsibilities r at all LW1 positions on sequence n
        // n runs over all sequences
#pragma omp parallel for reduction(+:llikelihood)
        for( size_t n = 0; n < seqs_.size(); n++ ){

            size_t 	L = seqs_[n]->getL();
            size_t 	LW1 = L - W_ + 1;
            size_t*	kmer = seqs_[n]->getKmer();
            float 	normFactor = 1.0f - q_;

            // initialize r_[n][i] and pos_[n][i]
            float pos_i = q_ / static_cast<float>( LW1 );
            for( size_t idx = 0; idx < ri[n].size(); idx++ ){
                r_[n][ri[n][idx]] = 1.0f;
                pos_[n][ri[n][idx]] = pos_i;
            }

            // update r based on the log odd ratios
            for( size_t idx = 0; idx < ri[n].size(); idx++ ){
                for( size_t j = 0; j < W_; j++ ){
                    size_t y = kmer[L-W_-ri[n][idx]+j] % Y_[K_+1];
                    r_[n][ri[n][idx]] *= s_[y][j];
                }
                // calculate the responsibilities and sum them up
                r_[n][ri[n][idx]] *= pos_[n][LW1-ri[n][idx]];
                normFactor += r_[n][ri[n][idx]];
            }

            // normalize responsibilities
            r_[n][0] /= normFactor;
            for( size_t idx = 0; idx < ri[n].size(); idx++ ){
                r_[n][ri[n][idx]] /= normFactor;
            }

            // for the unused positions
            for( size_t i = LW1; i < L; i++ ){
                r_[n][i] = 0.0f;
            }

            // calculate log likelihood over all sequences
            llikelihood += logf( normFactor );
        }

        llikelihood_ = llikelihood;
        /**
         * M-step for f_% motif occurrences
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
        // using the f_% r's
        // n runs over all sequences
        for( size_t n = 0; n < seqs_.size(); n++ ){
            size_t  L = seqs_[n]->getL();
            size_t* kmer = seqs_[n]->getKmer();

            for( size_t idx = 0; idx < ri[n].size(); idx++ ){
                for( size_t j = 0; j < W_; j++ ){
                    size_t y = kmer[L-W_-ri[n][idx]+j] % Y_[K_+1];
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
        if( verbose_ ) std::cout << iteration << "th iteration, delta_llikelihood=" << llikelihood_diff << std::endl;
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

    float N1 = 0.f;                   // expectation value of the count of sequences with at least a query motif

    for( size_t n = 0; n < seqs_.size(); n++ ){
        for( size_t i = 0; i < seqs_[n]->getL() - W_ + 1; i++ ){
            N1 += r_[n][i];
        }
    }

    q_ = ( seqs_.size() - N1 + 1.f ) / ( ( float )seqs_.size() + 2.f );

    if( verbose_ ) std::cout << "optimized q=" << q_ << std::endl;

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
        for( size_t y = 0; y < Y_[K_+1]; y++ ){
            std::cout << std::setprecision( 3 ) << n_[K_][y][j] << '\t';
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
        L = ss ? L : ( L - 1 ) / 2;

        for( size_t i = 0; i < seqs_[n]->getL()-W_+1; i++ ){

            if( r_[n][seqs_[n]->getL() -W_-i] >= cutoff ){
                ofile_pos << '>' << seqs_[n]->getHeader() << '\t' << L << '\t'
                          << ( ( i < L ) ? '+' : '-' ) << '\t' << i + 1 << ".." << i+W_ << '\t';
                for( size_t b = i; b < i+W_; b++ ){
                    ofile_pos << Alphabet::getBase( seqs_[n]->getSequence()[b] );
                }
                ofile_pos << std::endl;
            }
        }
    }

	// output responsibilities r[n][i]
	bool write_r = false;		// flag for writing out the responsibilities
	if( write_r ){
		std::string opath_r = opath + ".weights";
		std::ofstream ofile_r( opath_r.c_str() );
		for( size_t n = 0; n < seqs_.size(); n++ ){
			for( size_t i = 0; i < seqs_[n]->getL() - W_ + 1; i++ ){
				ofile_r << std::setprecision( 2 ) << r_[n][seqs_[n]->getL()-W_-i] << ' ';
			}
			ofile_r << std::endl;
		}
	}
}