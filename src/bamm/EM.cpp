#include "EM.h"

EM::EM( Motif* motif, BackgroundModel* bg, std::vector<int> folds ){

	motif_ = motif;
	bg_ = bg;

	int y, k, j, L;
	int W = motif_->getW();

	for( k = 0; k < std::max( Global::modelOrder+2, 4 ); k++ ){	// 4 is for cases when modelOrder < 2
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
		if( Global::setSlow ){
			L = posSeqs_[n]->getL() - W + 1;
		} else {
			L = posSeqs_[n]->getL();
		}
		r_[n] = ( float* )calloc( L, sizeof( float ) );
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

	long timestamp = time( NULL );

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
				if( Global::logEM )
					// calculate s in log space:
					s_[y][j] = logf( motif_->getV()[K_model][y][j] ) - logf( bg_->getV()[std::min( K_model, K_bg )][y_bg] );
				else
					// calculate s in linear space:
					s_[y][j] = motif_->getV()[K_model][y][j] / bg_->getV()[std::min( K_model, K_bg )][y_bg];
			}
		}

		// E-step: calculate posterior
		EStep();

		// M-step: update parameters
		MStep();

		// * optional: optimize parameter alpha
		if( !Global::noAlphaOptimization ){

		    // only run alpha optimization every 5th em-iteration.
		    if( EMIterations_ % Global::alphaIter == 0 ){
				testAlphaLearning();
				printf( "FINISHED!\n");
				exit(0);
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
			std::cout << EMIterations_ << "th iteration:	";
			std::cout << "para_diff = " << v_diff << ",	";
			std::cout << "log likelihood = " << llikelihood_ << " 	";
			if( llikelihood_diff < 0 && EMIterations_ > 1) std::cout << " decreasing...";
			std::cout << std::endl;
		}
		if( v_diff < Global::epsilon )					iterate = false;
		//if( llikelihood_diff < 0 && EMIterations_ > 1 )	iterate = false;

		// * testing: write out alpha, qfunc, gradient and posterior value for current EM iterations
		if( Global::TESTING ){
		    std::string opath = std::string( Global::outputDirectory ) + '/'
		            + std::string( Global::posSequenceBasename );
		    std::string opath_testing = opath + ".TESTING";
		    std::ofstream ofile_testing;
		    ofile_testing.open( opath_testing.c_str() , std::ios_base::app);

		    for( int k = 0; k < Global::modelOrder+1; k++ ){
		        ofile_testing << std::setprecision( 3 ) << alpha_[k][0] << ' ';
		    }
		    ofile_testing << std::scientific << calculateQfunc_gradient(alpha_[Global::modelOrder][0],Global::modelOrder) << ' ';
		    ofile_testing << std::scientific << calculateQfunc() << ' ';
		    ofile_testing << std::scientific << calculateLogPosterior() << ' ';
		    ofile_testing << std::endl;
		}
	}

	// calculate probabilities
	motif_->calculateP();

	// free memory
	for( y = 0; y < Y_[Global::modelOrder+1]; y++ ){
		free( v_prev[y] );
	}
	free( v_prev );

	fprintf( stdout, "\n--- Runtime for EM: %ld seconds ---\n", time( NULL )-timestamp );

    return 0;
}

void EM::EStep(){

	int L, LW1, k, y, j, i;
	float prior_i;											// positional preference prior state for p(z_n = i), motif present
	float prior_0 = 1 - q_;									// positional preference prior state for p(z_n = 0), no motif present
	float normFactor;										// normalize responsibilities r[n][i]
	int K_model = Global::modelOrder;
	int W = motif_->getW();
	llikelihood_ = 0.0f;									// reset log likelihood

	// calculate responsibilities r_[n][i] at position i in sequence n
	if( Global::setSlow ){
		// slow code:
		for( size_t n = 0; n < posSeqs_.size(); n++ ){		// n runs over all sequences
			L = posSeqs_[n]->getL();
			LW1 = L - W + 1;
			normFactor = 0.0f;								// reset normalization factor

			// reset r_[n][i]
			for( i = 0; i < LW1; i++ ){
				if( Global::logEM )
					// calculation in log space:
					r_[n][i] = 0.0f;
				else
					// calculation in linear space:
					r_[n][i] = 1.0f;
			}

			// when p(z_n > 0)
			prior_i = q_ / static_cast<float>( LW1 );		// p(z_n = i), i > 0
			for( i = 0; i < LW1; i++ ){
				for( j = 0; j < W; j++ ){
					k = std::min( i+j, K_model );
					y = posSeqs_[n]->extractKmer( i+j, k );
					if( y != -1 ){							// skip 'N' and other unknown alphabets
						if( Global::logEM )
							// calculation in log space:
							r_[n][i] += s_[y][j];
						else
							// calculation in linear space:
							r_[n][i] *= s_[y][j];
					} /*else if ( j < K_model ){			// for N exists and j < K occasions
						y = posSeqs_[n]->extractKmer( i+j, j );
						if( y != -1 ){
							r_[n][i] += s_[y][j];
						} else {
							r_[n][i] = 0.0f;
							break;
						}
					}*/ else {
						r_[n][i] = 0.0f;
						break;
					}
				}
				if( r_[n][i] != 0.0f ){
					if( Global::logEM )
						// calculation in exponential space:
						r_[n][i] = expf( r_[n][i] ) * prior_i;
					else
						// calculation in linear space:
						r_[n][i] = r_[n][i] * prior_i;
				}
				normFactor += r_[n][i];
			}

			// when p(z_n = 0)
			normFactor += prior_0;

			for( i = 0; i < LW1; i++ ){						// responsibility normalization
				r_[n][i] /= normFactor;
			}

			llikelihood_ += logf( normFactor );

		}
	} else {
		// fast code:
		for( size_t n = 0; n < posSeqs_.size(); n++ ){		// n runs over all sequences
			L = posSeqs_[n]->getL();
			LW1 = L - W + 1;
			normFactor = 0.0f;								// reset normalization factor

			// reset r_[n][i]
			for( i = 0; i < L; i++ ){
				if( Global::logEM )
					// calculation in log space:
					r_[n][i] = 0.0f;
				else
					// calculation in linear space:
					r_[n][i] = 1.0f;
			}

			// when p(z_n > 0)
			prior_i = q_ / static_cast<float>( LW1 );		// p(z_n = i), i > 0
			for( i = 0; i < L; i++ ){						// i runs over all nucleotides in sequence
				k = std::min( i, K_model );
				y = posSeqs_[n]->extractKmer( i, k );		// extract (k+1)-mer y from positions (i-k,...,i)
				for( j = 0; j < std::min( W, i+1 ); j++ ){	// j runs over all motif positions
					if( y != -1 ){							// skip 'N' and other unknown alphabets
						if( Global::logEM )
							// calculation in log space:
							r_[n][L-i+j-1] += s_[y][j];
						else
							// calculation in linear space:
							r_[n][L-i+j-1] *= s_[y][j];
					}/* else if( j < K_model ){
						y = posSeqs_[n]->extractKmer( i, K_model-j );
						if( y != -1 ){
							r_[n][L-i+j-1] += s_[y][j];
						} else {
							r_[n][L-i+j-1] = 0.0f;
							break;
						}
					}*/ else {
						r_[n][L-i+j-1] = 0.0f;
						break;
					}
				}
			}
			for( i = W-1; i < L; i++ ){
				if( r_[n][i] != 0.0f ){
					if( Global::logEM )
						// calculation in exponential space:
						r_[n][i] = expf( r_[n][i] ) * prior_i;
					else
						// calculation in linear space:
						r_[n][i] = r_[n][i] * prior_i;
				}
				normFactor += r_[n][i];
			}

			// when p(z_n = 0)
			normFactor += prior_0;

			// normalize responsibilities
			for( i = W-1; i < L; i++ ){
				r_[n][i] /= normFactor;
			}

			llikelihood_ += logf( normFactor );
		}
	}
}

void EM::MStep(){

	int L, LW1, k, y, y2, j, i;
	int W = motif_->getW();

	// reset the fractional counts n
	for( k = 0; k < Global::modelOrder+1; k++ ){
		for( y = 0; y < Y_[k+1]; y++ ){
			for( j = 0; j < W; j++ ){
				n_[k][y][j] = 0.0f;
			}
		}
	}

	// compute fractional occurrence counts for the highest order K
	if( Global::setSlow ){
		// slow code:
		for( size_t n = 0; n < posSeqs_.size(); n++ ){		// n runs over all sequences
			L = posSeqs_[n]->getL();
			LW1 = L - W + 1;
			for( i = 0; i < LW1; i++ ){						// i runs over all nucleotides in sequence
				for( j = 0; j < W; j++ ){					// j runs over all motif positions
					k = std::min( i+j, Global::modelOrder );
					y = posSeqs_[n]->extractKmer( i+j, k );
					if( y != -1 ){							// skip 'N' and other unknown alphabets
						n_[Global::modelOrder][y][j] += r_[n][i];
					}
				}
			}
		}
	} else {
		// fast code:
		for( size_t n = 0; n < posSeqs_.size(); n++ ){		// n runs over all sequences
			L = posSeqs_[n]->getL();
			LW1 = L - W + 1;
			for( i = 0; i < L; i++ ){
				k =  std::min( i, Global::modelOrder );
				y = posSeqs_[n]->extractKmer( i, k );
				for( j = 0; j < std::min( W, i+1 ); j++ ){
					if( y != -1 && ( i-j ) < LW1 ){			// skip 'N' and other unknown alphabets
						n_[Global::modelOrder][y][j] += r_[n][L-i+j-1];
					}
				}
			}
		}
	}

	// compute fractional occurrence counts from higher to lower order
	for( k = Global::modelOrder; k > 0; k-- ){				// k runs over all orders
		for( y = 0; y < Y_[k+1]; y++ ){
			y2 = y % Y_[k];									// cut off the first nucleotide in (k+1)-mer
			for( j = 0; j < W; j++ ){
				n_[k-1][y2][j] += n_[k][y][j];
			}
		}
	}

	// update model parameters v[k][y][j]
	motif_->updateV( n_, alpha_ );

}

//typedef int( EM::*EMMemFn)(double a);

void EM::optimizeAlphas( float min_brent, float max_brent, float tolerance ){
    //	if( Global::debugMode){
    //		fprintf( stderr, "\n ---------------------------------- \n");
    //		fprintf( stderr, "\n Verbose Output for Alpha Learning: \n");
    //		fprintf( stderr, "\n ---------------------------------- \n\n");
    //	}

    // for now: only optimize highest alpha
    int k = Global::modelOrder;
    //for( int k = 0; k < Global::modelOrder+1; k++ ){
        float optim_alpha = zbrent( *this, &EM::calculateQfunc_gradient, min_brent, max_brent, tolerance, k );
        //      if( Global::debugMode ){
        //          fprintf( stderr, "\n alpha_%d : \t %0.4f -> \t %0.4f \n ", k, alpha_[k][0], optim_alpha );
        //      }

        // only update in case a root is bracketed
        if( optim_alpha > 0 ){
            for( int j = 0; j < motif_->getW() ; j++ ){
                alpha_[k][j] = optim_alpha;
            }
            motif_->updateVbyK( n_, alpha_, k );
        }
    //}
}

void EM::testAlphaLearning( ){

	float alpha, alpha_min = 1, alpha_max = 5000;
	std::vector<float> gradient, qvalue, posterior;

	for(int k = 0; k <= Global::modelOrder; k++ ){

	    std::string opath = std::string( Global::outputDirectory ) + '/'
	            + std::string( Global::posSequenceBasename );

	    std::stringstream alphaIter;
	    alphaIter << Global::alphaIter;

        std::stringstream kstring;
        kstring << k;


	    std::string opath_n = opath + "_emIter_" + alphaIter.str() + "_Order_" + kstring.str() + ".AlphaTesting";
	    std::ofstream ofile_n( opath_n.c_str() );

	    for( alpha = alpha_min; alpha < alpha_max; alpha++ ){
	        // update alpha
	        for( int j = 0; j < motif_->getW(); j++ ){
	            alpha_[k][j] = alpha;
	        }
	        // adjust v_s to new alpha
	        motif_->updateVbyK( n_, alpha_ ,k );
	        // calculate and store Alpha_Gradient_Qfunc_LogPosterior
	        ofile_n << std::scientific << alpha << ' ';
	        ofile_n << std::scientific << calculateQfunc_gradient( alpha, k ) << ' ';
	        ofile_n << std::scientific << calculateQfunc( k ) << ' ';
	        ofile_n << std::scientific << calculateLogPosterior( k ) << ' ';
	        ofile_n << std::endl;
	    }
	    // reset v's for initial alpha
	    for( int j = 0; j < motif_->getW(); j++ ){
	        alpha_[k][j] = Global::modelAlpha[k];
	    }
	    motif_->updateVbyK( n_, alpha_ ,k );
	}
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
	        if( Global::setSlow ){
	            L = posSeqs_[n]->getL() - W + 1;
	        } else {
	            L = posSeqs_[n]->getL();
	        }
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
	                    if( Global::logEM )
	                        // calculation in log space:
	                        r_local[n][L-i+j-1] += s_[y][j];
	                    else
	                        // calculation in linear space:
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
	                if( Global::logEM )
	                    // calculation in exponential space:
	                    r_local[n][i] = expf( r_local[n][i] ) * prior_i;
	                else
	                    // calculation in linear space:
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
	float Qfunc = 0.0f;
	float prior_i, prior_0 = 1 - q_;

	int W = motif_->getW();
	int A = Alphabet::getSize();
	float*** v_motif = motif_->getV();
    float** v_bg = bg_->getV();

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
		Qfunc += ( r_[n][0] * logf( prior_0 ) + ( 1 - r_[n][0] ) * logf( prior_i ) );
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

	int k, y, j, i, L;
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

	// output responsibilities r[n][i]
	std::string opath_r = opath + "_emIter_" + alphaIter.str() + ".EMposterior";
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
