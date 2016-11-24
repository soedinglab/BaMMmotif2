#include "EM.h"

EM::EM( Motif* motif, BackgroundModel* bg, std::vector<int> folds ){

	motif_ = motif;
	bg_ = bg;

	int y, k, j, LW1;
	int W = motif_->getW();

	for( k = 0; k < std::max( Global::modelOrder+2,  Global::bgModelOrder+2 ); k++ ){
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
/*	for( size_t f_idx = 0; f_idx < folds_.size(); f_idx++ ){
		for( size_t s_idx = 0; s_idx < Global::posFoldIndices[folds_[f_idx]].size(); s_idx++ ){
			posSeqs_.push_back( Global::posSequenceSet->getSequences()[Global::posFoldIndices[folds_[f_idx]][s_idx]] );
		}
	}*/
	posSeqs_ = Global::posSequenceSet->getSequences();

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

	// allocate memory for motif position z_[n]
	z_ = ( int* )calloc( posSeqs_.size(), sizeof( int ) );
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

	free( z_ );

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
	float** v_before;										// hold the parameters of the highest-order before EM
	float q_func_old, q_func_new, l_post_old, l_post_new, l_prior_new, prev_alpha;


	// allocate memory for parameters v[y][j] with the highest order
	v_before = ( float** )calloc( Y_[Global::modelOrder+1], sizeof( float* ) );
	for( y = 0; y < Y_[Global::modelOrder+1]; y++ ){
		v_before[y] = ( float* )calloc( W, sizeof( float ) );
	}

    l_post_old = calculateLogPosterior();

	// iterate over
	while( iterate && ( EMIterations_ < Global::maxEMIterations ) ){

		EMIterations_++;

		// get parameter variables with highest order before EM
		for( y = 0; y < Y_[Global::modelOrder+1]; y++ ){
			for( j = 0; j < W; j++ ){
				v_before[y][j] = motif_->getV()[Global::modelOrder][y][j];
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

		// M-step: update model parameters
		MStep();

		// check EM-criteria
		q_func_old = calculateQfunc();

		// update model parameters v[k][y][j]
		if( Global::fixPseudos ){
			motif_->updateVbyK( n_, alpha_, K_model );
		} else {
			motif_->updateV( n_, alpha_ );
		}

		// * optional: optimize parameter alpha
    	prev_alpha = alpha_[K_model][0];
		if( !Global::noAlphaOptimization ){
		    // only run alpha optimization every xth em-iteration.
            if( EMIterations_ % Global::alphaIter == 0 && EMIterations_ > 10){
                optimizeAlphas();
//                std::cout << "optimAlpha = " << alpha_[Global::modelOrder][0] << "\n ";

            }
//		    if( EMIterations_ == 18 ){
//                optimizeAlphas();
//                std::cout << "LastoptimAlpha = " << alpha_[Global::modelOrder][0] << "\n ";
//                //testAlphaLearning();
//		        //std::cout << "LastoptimAlpha = " << alpha_[Global::modelOrder][0] << "\n ";
//
//		        exit(0);
//
//		    }
    		// check EM-criteria
    		q_func_new = calculateQfunc();
    		l_post_new = calculateLogPosterior();
    		l_prior_new = calculateLogPriors();

    		// check Qfunc increase
    		if( (q_func_new - q_func_old) < 0 ){
    			// reset alpha:
    			alpha_[K_model][0] = prev_alpha;
    			//reset V's:
    			motif_->updateVbyK( n_, alpha_, K_model );
    			// writeOut Qfunction Values
    			testAlphaLearning();
    		}
		}

		// * optional: optimize parameter q
		if( !Global::noQOptimization )		optimizeQ();

		// check parameter difference for convergence
		v_diff = 0.0f;
		for( y = 0; y < Y_[K_model+1]; y++ ){
			for( j = 0; j < W; j++ ){
				v_diff += fabsf( motif_->getV()[K_model][y][j] - v_before[y][j] );
			}
		}
		// check likelihood for convergence
		llikelihood_diff = llikelihood_ - llikelihood_prev;

		if( Global::verbose ){
			std::cout << EMIterations_ << " iteration:	";
			std::cout << "para_diff = " << v_diff << ",	";
			std::cout << "log likelihood = " << llikelihood_ << " 	";
			std::cout << "logPosterior " << l_post_new << "     ";
			if( llikelihood_diff < 0 && EMIterations_ > 1) std::cout << " decreasing... ";
			if( ( q_func_new - q_func_old ) < 0 ) std::cout << " ! qfunc decr.. !  ";
			if( ( l_post_new - l_post_old ) < 0 ) std::cout << " ! lPost decr.. !  ";
			std::cout << std::endl;
		}
		l_post_old = l_post_new;

		if( v_diff < Global::epsilon )					iterate = false;
		if( llikelihood_diff < 0 && EMIterations_ > 1 )	iterate = false;

		// * testing: write out alpha, qfunc, gradient and posterior value for current EM iterations
		if( Global::TESTING ){
		    std::stringstream alphaIter;
		    alphaIter << Global::alphaIter;

		    std::string opath = std::string( Global::outputDirectory ) + '/'
		            + std::string( Global::posSequenceBasename );
		    std::string opath_testing = opath + "emIter" + alphaIter.str() + ".TESTING";
		    std::ofstream ofile_testing;
		    ofile_testing.open( opath_testing.c_str() , std::ios_base::app);
		    ofile_testing << std::scientific << calculateQfunc_gradient(alpha_[Global::modelOrder][0],Global::modelOrder) << ' ';
		    ofile_testing << std::scientific << q_func_new - q_func_old << ' ';
		    ofile_testing << std::scientific << l_post_new << ' ';
		    ofile_testing << std::scientific << l_prior_new << ' ';
		    ofile_testing << std::scientific << llikelihood_ << ' ';

		    for( int k = 0; k < Global::modelOrder+1; k++ ){
		        ofile_testing << std::setprecision( 3 ) << alpha_[k][0] << ' ';
		    }
		    ofile_testing << std::endl;
		}
	}

	// calculate z[n]
	for( size_t n = 0; n < posSeqs_.size(); n++ ){
		int LW1 = posSeqs_[n]->getL() - W + 1;
		float max = 0;
		for( int i = 0; i < LW1; i++ ){
			if( r_[n][i] > max ){
				max = r_[n][i];
				z_[n] = LW1-i-1;
			}
		}
	}

	// reset n_z_[k][y][j]
	for( int k = 0; k < K_model+1; k++ ){
		for( y = 0; y < Y_[k+1]; y++ ){
			for( j = 0; j < W; j++ ){
				n_[k][y][j] = 0;
			}
		}
	}
	for( int i = 0; i < posSeqs_.size(); i++ ){
		for( j = 0; j < W; j++ ){
			y = posSeqs_[i]->extractKmer( z_[i]+j, std::min( z_[i]+j, K_model ) );
			n_[K_model][y][j]++;
		}
	}
	// calculate nz for lower order k
	for( int k = K_model; k > 0; k-- ){				// k runs over all lower orders
		for( y = 0; y < Y_[k+1]; y++ ){
			int y2 = y % Y_[k];					// cut off the first nucleotide in (k+1)-mer
			for( j = 0; j < W; j++ ){
				n_[k-1][y2][j] += n_[k][y][j];
			}
		}
	}
	// print kmer counts out
	std::cout << std::endl;
	for( j = 0; j < W; j++ ){
	for( int k = 0; k < K_model+1; k++ ){
		for( y = 0; y < Y_[k+1]; y++ ){

				std::cout << n_[k][y][j] << '\t';
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;


	// only for testing:
	std::cout << std::endl;
	for( size_t n = 0; n < 50/*posSeqs_.size()*/; n++ ){
		std::cout << z_[n] << '\t';
	}
	std::cout << std::endl;

	// calculate probabilities
	motif_->calculateP();

	// free memory
	for( y = 0; y < Y_[Global::modelOrder+1]; y++ ){
		free( v_before[y] );
	}
	free( v_before );

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
}

void EM::optimizeAlphas( float min_brent, float max_brent, float tolerance ){
   //for( int k = 1; k < Global::modelOrder+1; k++ ){
	int k = Global::modelOrder;
	float optim_alpha = zbrent( *this, &EM::calculateQfunc_gradient, min_brent, max_brent, tolerance, k );

	// only update in case a root is bracketed
	if( optim_alpha > 0 ){
		for( int j = 0; j < motif_->getW() ; j++ ){
			alpha_[k][j] = optim_alpha;
		}
		if( Global::fixPseudos ){
			motif_->updateVbyK( n_, alpha_, k );
		} else {
			motif_->updateV( n_, alpha_ );
		}
	}
}

void EM::testAlphaLearning( ){
    std::cout << "Starting AlphaLearningtesting...  ";
    std::cout << std::endl;
    float alpha, alpha_min = 1, alpha_max = 2e4;

    int k = Global::modelOrder;
        std::cout << k << " th Order " << ' ';
        std::cout << std::endl;
        std::string opath = std::string( Global::outputDirectory ) + '/'
                + std::string( Global::posSequenceBasename );
        std::stringstream emIter;
        emIter << EMIterations_;
        std::stringstream kstring;
        kstring << k;
        std::string opath_n = opath + "_emIter_" + emIter.str() + "_Order_" + kstring.str() + ".AlphaTesting";
        std::ofstream ofile_n( opath_n.c_str() );

        for( alpha = alpha_min; alpha < alpha_max; alpha++ ){
            std::cout << "  alpha= " << alpha << ' ';
            // update alpha
            for( int j = 0; j < motif_->getW(); j++ ){
                alpha_[k][j] = alpha;
            }
            // adjust v_s to new alpha
            motif_->updateVbyK( n_, alpha_ ,k );
            // calculate and store Alpha_Gradient_Qfunc_LogPosterior
            ofile_n << std::scientific << alpha << ' ';
            ofile_n << std::scientific << calculateQfunc_gradient( alpha, k ) << ' ';
            std::cout << "  gradient= " << calculateQfunc_gradient( alpha, k ) << ' ';
            ofile_n << std::scientific << calculateQfunc( k ) << ' ';
            std::cout << "  Qfunction= " << calculateQfunc( k ) ;
            ofile_n << std::scientific << calculateLogPosterior( k ) << ' ';
            ofile_n << std::endl;
            std::cout << std::endl;
        }
        std::cout << "              Resetting v's " ;
        std::cout << std::endl << std::flush;
        // reset v's for initial alpha
        for( int j = 0; j < motif_->getW(); j++ ){
            alpha_[k][j] = Global::modelAlpha[k];
        }
        motif_->updateVbyK( n_, alpha_ ,k );
}

void EM::optimizeQ(){

	// optimize hyper-parameter q
	// motif.updateV()
}

float EM::calculateLogPriors( int K ){

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
			y2 = y % Y_[K];                         // cut off the first nucleotide in (k+1)-mer y
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

float EM::calculateLogPosterior( int K ){

	return llikelihood_ + calculateLogPriors( K );
}

float EM::calculateQfunc( int K ){

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
	Qfunc += calculateLogPriors( K );

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

// only for obtaining z for CGS testing
int* EM::getZ(){
	return z_;
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
				ofile_n << ( int )n_[k][y][j] << ' ';
			}
			ofile_n << std::endl;
		}
		ofile_n << std::endl;
	}

/*
	//TODO: this file is too large for benchmarking
	// output responsibilities r[n][i]
	std::string opath_r = opath + ".EMposterior";
	std::ofstream ofile_r( opath_r.c_str() );
	for( size_t n = 0; n < posSeqs_.size(); n++ ){
		int L = posSeqs_[n]->getL();
		for( int i = L-1; i > W-2; i-- ){
			ofile_r << std::scientific << r_[n][i] << ' ';
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
