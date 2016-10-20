#include "EM.h"

EM::EM( Motif* motif, BackgroundModel* bg, std::vector<int> folds ){

	motif_ = motif;
	bg_ = bg;

	int y, k, j, LW1;
	int W = motif_->getW();

	for( k = 0; k < std::max( Global::modelOrder+2,  Global::bgModelOrder+2 ); k++ ){	// 4 is for cases when modelOrder < 2
		Y_.push_back( ipow( Alphabet::getSize(), k ) );
	}

	if( folds.empty() ){										// for the complete sequence set
		folds_.resize( Global::cvFold );
		std::iota( std::begin( folds_ ), std::end( folds_ ), 0 );
	} else {													// for cross-validation
		folds_ = folds;
	}

	// copy selected positive sequences due to folds
	posSeqs_.clear();
	for( size_t f_idx = 0; f_idx < folds_.size(); f_idx++ ){
		for( size_t s_idx = 0; s_idx < Global::posFoldIndices[folds_[f_idx]].size(); s_idx++ ){
			posSeqs_.push_back( Global::posSequenceSet->getSequences()[Global::posFoldIndices[folds_[f_idx]][s_idx]] );
		}
	}

	// allocate memory for s_[y][j] and initialize it
	s_ = ( float** )calloc( Y_[Global::modelOrder+1], sizeof( float* ) );
	for( y = 0; y < Y_[Global::modelOrder+1]; y++ ){
		s_[y] = ( float* )calloc( W, sizeof( float ) );
	}

	// allocate memory for r_[n][i] and initialize it
	r_ = ( float** )calloc( posSeqs_.size(), sizeof( float* ) );
	for( size_t n = 0; n < posSeqs_.size(); n++ ){
		LW1 = posSeqs_[n]->getL() - W + 1;
		r_[n] = ( float* )calloc( LW1, sizeof( float ) );
	}

	// allocate memory for n_[k][y][j] and probs_[k][y][j] and initialize them
	n_ = ( float*** )calloc( Global::modelOrder+1, sizeof( float** ) );
	for( k = 0; k < Global::modelOrder+1; k++ ){
		n_[k] = ( float** )calloc( Y_[k+1], sizeof( float* ) );
		for( y = 0; y < Y_[k+1]; y++ ){
			n_[k][y] = ( float* )calloc( W, sizeof( float ) );
		}
	}

	// allocate memory for alpha_[k][j] and initialize it
	alpha_ = ( float** )malloc( ( Global::modelOrder+1 ) * sizeof( float* ) );
	for( k = 0; k < Global::modelOrder+1; k++ ){
		alpha_[k] = ( float* )malloc( W * sizeof( float ) );
		for( j = 0; j < W; j++ ){
			alpha_[k][j] = Global::modelAlpha[k];
		}
	}
}

EM::~EM(){

	for( int y = 0; y < Y_[Global::modelOrder+1]; y++ ){
		free( s_[y] );
	}
	free( s_ );

	for( size_t n = 0; n < posSeqs_.size(); n++ ){
		free( r_[n] );
	}
	free( r_ );

	for( int k = 0; k < Global::modelOrder+1; k++ ){
		for( int y = 0; y < Y_[k+1]; y++ ){
			free( n_[k][y] );
		}
		free( n_[k] );
	}
	free( n_ );

	for( int k = 0; k < Global::modelOrder+1; k++ ){
		free( alpha_[k] );
	}
	free( alpha_ );

}
void EM::testFunctions(){

    // test if Qfunction fits to its gradient
    float change = 0.001f;
    for(int i = 0; i < 5; i++){
        // 1. calculate Qfunction and Gradient for start alpha
        float Q_old = calculateQfunc();
        float Grad_old = calculateQfunc_gradient( alpha_[Global::modelOrder][0], Global::modelOrder );

        // 2. change alpha parameter slightly
        for(int j = 0; j < motif_->getW(); j++ ){
            alpha_[Global::modelOrder][j] += change;
        }

        // 3. calculate Qfunction and Gradient with new alphas
        float Q_new = calculateQfunc();
        float Grad_new = calculateQfunc_gradient( alpha_[Global::modelOrder][0], Global::modelOrder );

        // 4. estimate Gradient approximation
        float Grad_approx =  ( Q_new - Q_old )/ change;

        fprintf(stdout, "Change \t Q_old \t Q_new \t Grad_Aprox \t Grad_old \t Grad_new \n");
        fprintf(stdout, "%f \t %f \t %f \t %f \t %f \t %f \n\n", change, Q_old, Q_new, Grad_approx, Grad_old, Grad_new );

        // 5. restore original alphas
        for(int j = 0; j < motif_->getW(); j++ ){
            alpha_[Global::modelOrder][j] -= change;
        }

        change = change*10;
    }
}

int EM::learnMotif(){

	printf( " ______\n"
			"|      |\n"
			"|  EM  |\n"
			"|______|\n\n" );

	clock_t t0 = clock();
	bool iterate = true;									// flag for iterating before convergence
	int W = motif_->getW();
	int K_model = Global::modelOrder;
	int K_bg = Global::bgModelOrder;

	int y, y_bg, j;
	float v_diff, llikelihood_prev, llikelihood_diff = 0.0f;
	float** v_prev;											// hold the parameters of the highest-order before EM

	// allocate memory for parameters v[y][j] with the highest order
	v_prev = ( float** )calloc( Y_[Global::modelOrder+1], sizeof( float* ) );
	for( y = 0; y < Y_[Global::modelOrder+1]; y++ ){
		v_prev[y] = ( float* )calloc( W, sizeof( float ) );
	}



	// iterate over
	while( iterate && ( EMIterations_ < Global::maxEMIterations ) ){

		EMIterations_++;

		// before EM, get parameter variables before EM
		for( y = 0; y < Y_[Global::modelOrder+1]; y++ ){
			for( j = 0; j < W; j++ ){
				v_prev[y][j] = motif_->getV()[Global::modelOrder][y][j];
			}
		}
		llikelihood_prev = llikelihood_;

		// compute log odd scores s[y][j], log likelihoods of the highest order K
		for( y = 0; y < Y_[K_model+1]; y++ ){
			for( j = 0; j < W; j++ ){
				y_bg = y % Y_[K_bg+1];
				s_[y][j] = motif_->getV()[K_model][y][j] / bg_->getV()[std::min( K_model, K_bg )][y_bg];
			}
		}

		// E-step: calculate posterior
		EStep();

		// M-step: update parameters
		MStep();

		// * optional: optimize parameter alpha
		if( !Global::noAlphaOptimization ){
		    // only run alpha optimization every xth em-iteration.
		    if( EMIterations_ % Global::alphaIter == 0 && EMIterations_ > 10){
		        optimizeAlphas();
		    }
		}

		// * optional: optimize parameter q
		if( !Global::noQOptimization )		optimizeQ();

				// check parameter difference for convergence
		v_diff = 0.0f;

		for( y = 0; y < Y_[Global::modelOrder+1]; y++ ){
			for( j = 0; j < W; j++ ){
				v_diff += fabsf( motif_->getV()[Global::modelOrder][y][j] - v_prev[y][j] );
			}
		}
		// check likelihood for convergence
		llikelihood_diff = llikelihood_ - llikelihood_prev;

		if( Global::verbose ){
			std::cout << EMIterations_ << " iteration:	";
			std::cout << "para_diff = " << v_diff << ",	";
			std::cout << "log likelihood = " << llikelihood_ << " 	";
			if( llikelihood_diff < 0 && EMIterations_ > 1) std::cout << " decreasing... ";
			std::cout << std::endl;
		}

		if( v_diff < Global::epsilon )					iterate = false;
		if( llikelihood_diff < 0 && EMIterations_ > 1 )	iterate = false;

		// * testing: write out alpha, qfunc, gradient and posterior value for current EM iterations
	}

	// calculate probabilities
	motif_->calculateP();

	// free memory
	for( y = 0; y < Y_[Global::modelOrder+1]; y++ ){
		free( v_prev[y] );
	}
	free( v_prev );

	fprintf( stdout, "\n--- Runtime for EM: %.4f seconds ---\n", ( ( float )( clock() - t0 ) ) / CLOCKS_PER_SEC );

    return 0;
}

void EM::EStep(){

	float prior_i;												// positional preference prior state for p(z_n = i), motif present
	float prior_0 = 1 - q_;										// positional preference prior state for p(z_n = 0), no motif present
	int K_model = Global::modelOrder;
	int W = motif_->getW();
	llikelihood_ = 0.0f;										// reset log likelihood

	// calculate responsibilities r_[n][i] at position i in sequence n
	for( size_t n = 0; n < posSeqs_.size(); n++ ){				// n runs over all sequences
		int L = posSeqs_[n]->getL();
		int LW1 = L - W + 1;
		float normFactor = 0.0f;								// reset normalization factor

		// reset r_[n][i]
		for( int i = 0; i < LW1; i++ )	r_[n][i] = 1.0f;

		// when p(z_n > 0)
		prior_i = q_ / static_cast<float>( LW1 );				// p(z_n = i), i > 0
		for( int ij = 0; ij < L; ij++ ){						// ij = i+j runs over all positions in sequence
			int y = posSeqs_[n]->extractKmer( ij, std::min( ij, K_model ) );	// extract (k+1)-mer y from positions (i-k,...,i)
			for( int j = std::max( 0, ij-L+W ); j < std::min( W, ij+1 ); j++ ){	// j runs over all motif positions
				if( y != -1 ){									// skip 'N' and other unknown alphabets
					r_[n][L-W-ij+j] *= s_[y][j];
				} else {
					r_[n][L-W-ij+j] = 0.0f;
					break;
				}
			}
		}

		for( int i = 0; i < LW1; i++ ){
			if( r_[n][i] != 0.0f ){
				r_[n][i] *= prior_i;
			}
			normFactor += r_[n][i];
		}

		// when p(z_n = 0), r = 1 - q_
		normFactor += prior_0;

		// normalize responsibilities
		for( int i = 0; i < LW1; i++ ){
			r_[n][i] /= normFactor;
		}

		llikelihood_ += logf( normFactor );
	}
}

void EM::MStep(){

	int W = motif_->getW();

	// reset the fractional counts n
	for( int k = 0; k < Global::modelOrder+1; k++ ){
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
			int y = posSeqs_[n]->extractKmer( ij, std::min( ij, Global::modelOrder ) );
			for( int j = std::max( 0, ij-L+W ); j < std::min( W, ij+1 ); j++ ){
				if( y != -1 && ( ij-j ) < LW1 ){				// skip 'N' and other unknown alphabets
					n_[Global::modelOrder][y][j] += r_[n][L-W-ij+j];
				}
			}
		}
	}

	// compute fractional occurrence counts from higher to lower order
	for( int k = Global::modelOrder; k > 0; k-- ){				// k runs over all orders
		for( int y = 0; y < Y_[k+1]; y++ ){
			int y2 = y % Y_[k];									// cut off the first nucleotide in (k+1)-mer
			for( int j = 0; j < W; j++ ){
				n_[k-1][y2][j] += n_[k][y][j];
			}
		}
	}

	// update model parameters v[k][y][j]
	motif_->updateV( n_, alpha_ );

}

//typedef int( EM::*EMMemFn)(double a);

void EM::optimizeAlphas( float min_brent, float max_brent, float tolerance ){
       for( int k = 1; k < Global::modelOrder+1; k++ ){
        float optim_alpha = zbrent( *this, &EM::calculateQfunc_gradient, min_brent, max_brent, tolerance, k );

        // only update in case a root is bracketed
        if( optim_alpha > 0 ){
            for( int j = 0; j < motif_->getW() ; j++ ){
                alpha_[k][j] = optim_alpha;
            }
            motif_->updateV( n_, alpha_ );
        }
    }
}

void EM::testAlphaLearning( ){
}

void EM::optimizeQ(){

	// optimize hyper-parameter q
	// motif.updateV()
}

float EM::calculateLogPosterior( int K ){
	int L, LW1, i,j,k,y,y2;
	float prior_i;
	float prior_0 = 1-q_;
	float lPosterior = 0.0f;
	float normFactor = 0.0f;
    int A = Alphabet::getSize();
	int W = motif_->getW();
	float*** v_motif = motif_->getV();
	float** v_bg = bg_->getV();


	float** r_local = ( float** )calloc( posSeqs_.size(), sizeof( float* ) );
	    for( size_t n = 0; n < posSeqs_.size(); n++ ){
	        L = posSeqs_[n]->getL();
	        r_local[n] = ( float* )calloc( L, sizeof( float ) );
	    }

	    for( size_t n = 0; n < posSeqs_.size(); n++ ){      // n runs over all sequences
	        L = posSeqs_[n]->getL();
	        LW1 = L - W + 1;
	        normFactor = 0.0f;                              // reset normalization factor

	        // when p(z_n > 0)
	        prior_i = q_ / static_cast<float>( LW1 );       // p(z_n = i), i > 0
	        for( i = 0; i < L; i++ ){                       // i runs over all nucleotides in sequence
	            k = std::min( i, K );
	            y = posSeqs_[n]->extractKmer( i, k );       // extract (k+1)-mer y from positions (i-k,...,i)
	            for( j = 0; j < std::min( W, i+1 ); j++ ){  // j runs over all motif positions
	                if( y != -1 ){                          // skip 'N' and other unknown alphabets
	                    r_local[n][L-i+j-1] *= s_[y][j];
	                }
	                else {
	                    r_local[n][L-i+j-1] = 0.0f;
	                    break;
	                }
	            }
	        }
	        for( i = W-1; i < L; i++ ){
	            if( r_local[n][i] != 0.0f ){
	                r_local[n][i] = r_local[n][i] * prior_i;
	            }
	            normFactor += r_local[n][i];
	        }
	        // when p(z_n = 0)
	        normFactor += prior_0;

	        lPosterior += logf(normFactor);
	    }

	    // the second and third parts of log Posterior Probability
	        for( j = 0; j < W; j++ ){
	            // the second part
	            lPosterior += ( float )Y_[K] * lgammaf( alpha_[K][j] + ( float )A );
	            // the second and third terms
	            for( y = 0; y < Y_[K+1]; y++ ){
	                // the second term
	                y2 = y % Y_[K];                         // cut off the first nucleotide in (k+1)-mer y
	                if( K == 0 ){
	                    lPosterior -= lgammaf( alpha_[K][j] * v_bg[K][y] + 1.0f );
	                    // the third term
	                    lPosterior += alpha_[K][j] * v_bg[K][y] * logf( v_motif[K][y][j] );
	                }
	                if( K > 0 ){
	                    lPosterior -= lgammaf( alpha_[K][j] * v_motif[K-1][y2][j] + 1.0f );
	                    // the third term
	                    lPosterior += alpha_[K][j] * v_motif[K-1][y2][j] * logf( v_motif[K][y][j] );
	                }
	            }
	            // the forth part
	            lPosterior += ( - 2.0f * logf( alpha_[K][j] ) - Global::modelBeta * powf( Global::modelGamma, ( float )K ) /
	                    alpha_[K][j] + logf( Global::modelBeta * powf( Global::modelGamma, ( float )K ) ) );
	        }

	return lPosterior;
}

float EM::calculateQfunc( int K ){

	int L, y, y2, j;
	float prior_i, prior_0 = 1 - q_;

	int W = motif_->getW();
	int A = Alphabet::getSize();
	float*** v_motif = motif_->getV();
    float** v_bg = bg_->getV();
	float Qfunc = 0.0f;
    // the first part of Q function
	for( y = 0; y < Y_[K+1]; y++ ){
		for( j = 0; j < W; j++ ){
			Qfunc += n_[K][y][j] * logf( v_motif[K][y][j] );
		}
	}

	// the second part of Q function
	for( size_t n = 0; n < posSeqs_.size(); n++ ){
		L = posSeqs_[n]->getL();
		prior_i = q_ / static_cast<float>( L - W + 1 );
		Qfunc += ( prior_0 * logf( prior_0 ) + q_ * logf( prior_i ) );
	}

	// the third and forth parts of Q function
	for( j = 0; j < W; j++ ){
		// the third part of Q function
		// the first term
		Qfunc += ( float )Y_[K] * lgammaf( alpha_[K][j] + ( float )A );
		// the second and third terms
		for( y = 0; y < Y_[K+1]; y++ ){
			// the second term
		    y2 = y % Y_[K];							// cut off the first nucleotide in (k+1)-mer y
		    if( K == 0 ){
		        Qfunc -= lgammaf( alpha_[K][j] * v_bg[K][y] + 1.0f );
		        // the third term
		        Qfunc += alpha_[K][j] * v_bg[K][y] * logf( v_motif[K][y][j] );
		    }
		    if( K > 0 ){
		        Qfunc -= lgammaf( alpha_[K][j] * v_motif[K-1][y2][j] + 1.0f );
		        // the third term
		        Qfunc += alpha_[K][j] * v_motif[K-1][y2][j] * logf( v_motif[K][y][j] );
		    }
		}

		// the forth part of Q function
		Qfunc += ( - 2.0f * logf( alpha_[K][j] ) - Global::modelBeta * powf( Global::modelGamma, ( float )K ) /
				alpha_[K][j] + logf( Global::modelBeta * powf( Global::modelGamma, ( float )K ) ) );
	}

	return Qfunc;
}

float EM::calculateQfunc_gradient( float alpha ,int K ){

	int j, y, y2;
	float Qfunc_grad = 0.0f;
    float*** v_motif = motif_->getV();
	float** v_bg = bg_->getV();
	int W = motif_->getW();
    int A = Alphabet::getSize();

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

void EM::print(){

	std::cout << " ____________________________________" << std::endl;
	std::cout << "|                                    |" << std::endl;
	std::cout << "| Conditional probabilities after EM |" << std::endl;
	std::cout << "|____________________________________|" << std::endl;
	int W = motif_->getW();
	for( int j = 0; j < W; j++ ){
		for( int k = 0; k < Global::modelOrder+1; k++ ){
			for( int y = 0; y < Y_[k+1]; y++ ){
				std::cout << std::scientific << motif_->getV()[k][y][j] << '\t';
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
}

void EM::write(){

	/**
	 * save EM parameters in four flat files:
	 * (1) posSequenceBasename.EMcounts:	refined fractional counts of (k+1)-mers
	 * (2) posSequenceBasename.EMposterior: responsibilities, posterior distributions
	 * (3) posSequenceBasename.EMalpha:		hyper-parameter alphas
	 * (4) posSequenceBasename.EMlogScores:	log scores
	 */

	int k, y, j;
	int W = motif_->getW();

	std::stringstream alphaIter;
	alphaIter << Global::alphaIter;

	std::string opath = std::string( Global::outputDirectory ) + '/'
						+ std::string( Global::posSequenceBasename );

	// output (k+1)-mer counts n[k][y][j]
	std::string opath_n = opath + "_emIter_" + alphaIter.str() + ".EMcounts";
	std::ofstream ofile_n( opath_n.c_str() );
	for( j = 0; j < W; j++ ){
		for( k = 0; k < Global::modelOrder+1; k++ ){
			for( y = 0; y < Y_[k+1]; y++ ){
				ofile_n << std::scientific << n_[k][y][j] << ' ';
			}
			ofile_n << std::endl;
		}
		ofile_n << std::endl;
	}

/*	TODO: this file is too large for benchmarking
	// output responsibilities r[n][i]
	int i, L;
	std::string opath_r = opath + ".EMposterior";
	std::ofstream ofile_r( opath_r.c_str() );
	for( size_t n = 0; n < posSeqs_.size(); n++ ){
		L = posSeqs_[n]->getL();
		if( Global::setSlow ){
			for( i = 0; i < L-W+1; i++ ){
				ofile_r << std::scientific << r_[n][i] << ' ';
			}
		} else {
			for( i = L-1; i > W-2; i-- ){
				ofile_r << std::scientific << r_[n][i] << ' ';
			}
		}
		ofile_r << std::endl;
	}
*/

	// output parameter alphas alpha[k][j]
	std::string opath_alpha = opath + "_emIter_" + alphaIter.str() + ".EMalpha";
	std::ofstream ofile_alpha( opath_alpha.c_str() );
	for( k = 0; k < Global::modelOrder+1; k++ ){
		ofile_alpha << "k = " << k << std::endl;
		for( j = 0; j < W; j++ ){
			ofile_alpha << std::setprecision( 3 ) << alpha_[k][j] << ' ';
		}
		ofile_alpha << std::endl;
	}

	// output log scores s[y][j]
	std::string opath_s = opath + "_emIter_" + alphaIter.str() + ".EMlogScores";
	std::ofstream ofile_s( opath_s.c_str() );
	for( y = 0; y < Y_[Global::modelOrder+1]; y++ ){
		for( j = 0; j < W; j++ ){
			ofile_s << std::fixed << std::setprecision( 6 ) << s_[y][j] << '	';
		}
		ofile_s << std::endl;
	}
}
