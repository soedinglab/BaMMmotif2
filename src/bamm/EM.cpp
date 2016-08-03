/*
 * EM.cpp
 *
 *  Created on: Apr 18, 2016
 *      Author: administrator
 */

#include "EM.h"

EM::EM( Motif* motif, BackgroundModel* bg, std::vector<int> folds ){

	motif_ = motif;
	v_motif_ = motif_->getV();
	W_ = motif_->getW();

	bg_ = bg;
	v_bg_ = bg_->getVbg();

	if( folds.size() != 0 )												// for cross-validation
		folds_ = folds;

	int n, y, k, j;
	// allocate memory for s_[y][j] and initialize it
	s_ = ( float** )calloc( Global::powA[K_+1], sizeof( float* ) );
	for( y = 0; y < Global::powA[K_+1]; y++ )
		s_[y] = ( float* )calloc( W_, sizeof( float ) );

    // allocate memory for r_[n][i] and initialize it
	r_ = ( float** )calloc( posSetN_, sizeof( float* ) );
    for( n = 0; n < posSetN_; n++ ){
    	if( Global::setSlow )
    		r_[n] = ( float* )calloc( posSeqs_[n].getL()-W_+2, sizeof( float ) );
    	else
    		r_[n] = ( float* )calloc( posSeqs_[n].getL()+1, sizeof( float ) );
    }
    // allocate memory for n_[k][y][j] and p_[k][y][j] and initialize them
	n_ = ( float*** )calloc( K_+1, sizeof( float** ) );
	p_ = ( float*** )calloc( K_+1, sizeof( float** ) );
	for( k = 0; k < K_+1; k++ ){
		n_[k] = ( float** )calloc( Global::powA[k+1], sizeof( float* ) );
		p_[k] = ( float** )calloc( Global::powA[k+1], sizeof( float* ) );
		for( y = 0; y < Global::powA[k+1]; y++ ){
			n_[k][y] = ( float* )calloc( W_, sizeof( float ) );
			p_[k][y] = ( float* )calloc( W_, sizeof( float ) );
		}
	}

	// allocate memory for alpha_[k][j] and initialize it
	alpha_ = ( float** )malloc( ( K_+1 ) * sizeof( float* ) );
	for( k = 0; k < K_+1; k++ ){
		alpha_[k] = ( float* )malloc( W_ * sizeof( float ) );
		for( j = 0; j < W_; j++ )										// initialize alpha_ with global alpha
			if( k == 0 )
				alpha_[k][j] = 1;
			else
				alpha_[k][j] = 20 * Global::ipow( 3, k-1 );
	}
}

EM::~EM(){

	int k, y, n;
	for( y = 0; y < Global::powA[K_+1]; y++ )
		free( s_[y] );
	free( s_ );

	for( n = 0; n < posSetN_; n++ )
		free( r_[n] );
	free( r_ );

	for( k = 0; k < K_+1; k++ ){
		for( y = 0; y < Global::powA[k+1]; y++ ){
			free( n_[k][y] );
			free( p_[k][y] );
		}
		free( n_[k] );
		free( p_[k] );
	}
	free( n_ );
	free( p_ );

	for( k = 0; k < K_+1; k++ )
		free( alpha_[k] );
	free( alpha_ );

	std::cout << "Destructor for EM class works fine. \n";

}

int EM::learnMotif(){

	printf( " ______\n"
			"|      |\n"
			"|  EM  |\n"
			"|______|\n\n" );

	bool iterate = true;												// flag for iterating before convergence
	int k, y, j, y2, yk;
	float v_diff, llikelihood_prior, llikelihood_diff, Qfunc_prior, Qfunc_diff;
	float** v_prior;													// hold the parameters before EM


	// allocate memory for prior and posterior parameters v[y][j] with the highest order
	v_prior = ( float** )calloc( Global::powA[K_+1], sizeof( float* ) );
	for( y = 0; y < Global::powA[K_+1]; y++ )
		v_prior[y] = ( float* )calloc( W_, sizeof( float ) );

	// iterate over
	while( iterate && ( EMIterations_ < Global::maxEMIterations ) ){

		EMIterations_++;

		// before EM, get the prior parameter variables
		for( y = 0; y < Global::powA[K_+1]; y++ )
			for( j = 0; j < W_; j++ )
				v_prior[y][j] = v_motif_[K_][y][j];
		llikelihood_prior = llikelihood_;
		Qfunc_prior = Qfunc_;

		// E-step: calculate posterior
		EStep();

		// M-step: update parameters
		MStep();

		// * optional: optimizeAlphas()
		if( !Global::noAlphaOptimization )	optimizeAlphas();

		// * optional: optimizeQ()
		if( !Global::noQOptimization )		optimizeQ();

		// * check likelihood for convergence
		v_diff = 0.0f;													// reset difference of posterior probabilities
		for( y = 0; y < Global::powA[K_+1]; y++ )
			for( j = 0; j < W_; j++ )
				v_diff += fabsf( v_motif_[K_][y][j] - v_prior[y][j] );
		llikelihood_diff = llikelihood_ - llikelihood_prior;
		Qfunc_diff = Qfunc_ - Qfunc_prior;

		if( Global::verbose ){
			std::cout << "At " << EMIterations_ << "th iteration:\n";
			std::cout << "	- parameter difference = " << v_diff << "\n";
			std::cout << "	- log likelihood = " << llikelihood_ << "\n";
			std::cout << "	- log likelihood difference = " << llikelihood_diff;
			if( llikelihood_diff < 0 ) std::cout << " decreasing...";
			std::cout << std::endl;
			std::cout << "	- Q function = " << Qfunc_ << "\n";
			std::cout << "	- Q function difference = " << Qfunc_diff;
			if( Qfunc_diff < 0 ) std::cout << " decreasing...";
			std::cout << std::endl;
		}

		if( v_diff < Global::epsilon )	iterate = false;
	}

/*
    // * remark: iterate over sequences using Global::posFoldIndices and folds_
	if( folds_.size() > 0 )
		for( unsigned int f = 0; f < folds_.size() ; f++ ){
			int fold = folds_[f];
			for( unsigned int n = 0; n < Global::posFoldIndices[fold].size(); n++ )
				Sequence seq = Global::posSequenceSet->getSequences()[Global::posFoldIndices[fold][n]];
		}
*/

//	if( Global::verbose ) print();

	// for calculating probabilities, i.e. p(ACG) = p(G|AC) * p(AC)
	// when k = 0:
	for( j = 0; j < W_; j++ )
		for( y = 0; y < Global::powA[1]; y++ )
			p_[0][y][j] = v_motif_[0][y][j];

	// when k > 0:
	for( k = 1; k < K_+1; k++){
		for( y = 0; y < Global::powA[k+1]; y++ ){
			y2 = y % Global::powA[k];									// cut off the first nucleotide in (k+1)-mer
			yk = y / Global::powA[1];									// cut off the last nucleotide in (k+1)-mer
			for( j = 0; j < k; j++ )
				p_[k][y][j] = p_[k-1][y2][j];							// i.e. p(ACG) = p(CG)
			for( j = k; j < W_; j++ )
				p_[k][y][j] =  v_motif_[k][y][j] * p_[k-1][yk][j-1];
		}
	}

	// free memory
	for( y = 0; y < Global::powA[K_+1]; y++ )
		free( v_prior[y] );
	free( v_prior );

    return 0;
}

void EM::EStep(){

	int L, LW1, y, j, y3, n, i;
	float prior_i;														// positional preference prior state for p(z_n = i), motif present
	float prior_0 = 1 - q_;												// positional preference prior state for p(z_n = 0), no motif present
	float normFactor;													// normalize responsibilities r[n][i]
	llikelihood_ = 0.0f;												// reset log likelihood

	// compute log odd scores s[y][j], log likelihoods with the highest order K
	for( y = 0; y < Global::powA[K_+1]; y++ ){
		for( j = 0; j < W_; j++ ){
			if( K_ <= 2 ){
				s_[y][j] = logf( v_motif_[K_][y][j] ) - logf( v_bg_[K_][y] );
			} else {
				y3 = y % Global::powA[3];								// 3 rightmost nucleotides in (k+1)-mer y
				s_[y][j] = logf( v_motif_[K_][y][j] ) - logf( v_bg_[k_bg_][y3] );
			}
		}
	}

	// calculate responsibilities r_[n][i] at position i in sequence n
	if( Global::setSlow ){
		// slow code:
		for( n = 0; n < posSetN_; n++ ){								// n runs over all sequences
			L = posSeqs_[n].getL();
			LW1 = L - W_ + 1;
			normFactor = 0.0f;											// reset normalization factor

			// reset r_[n][i]
			for( i = 1; i <= LW1; i++ )
				r_[n][i] = 0.0f;

			// when p(z_n > 0)
			prior_i = q_ / ( float )LW1;								// p(z_n = i), i > 0
			for( i = 1; i <= LW1; i++ ){
				for( j = 0; j < W_; j++ ){
					y = posSeqs_[n].extractKmer( i+j-1, ( j < K_ ) ? j : K_ );
					if( y != -1 )										// skip 'N' and other unknown alphabets
						r_[n][i] += s_[y][j];
					else
						break;
				}
				r_[n][i] = expf( r_[n][i] ) * prior_i;
				normFactor += r_[n][i];
			}

			// when p(z_n = 0)
			r_[n][0] = prior_0;
			normFactor += r_[n][0];

			for( i = 0; i <= LW1; i++ )								// responsibility normalization
				r_[n][i] /= normFactor;
			llikelihood_ += logf( normFactor );
		}
	} else {
		// fast code:
		for( n = 0; n < posSetN_; n++ ){								// n runs over all sequences
			L = posSeqs_[n].getL();
			LW1 = L - W_ + 1;
			normFactor = 0.0f;											// reset normalization factor

			// reset r_[n][i]
			for( i = 0; i <= L; i++ )
				r_[n][i] = 0.0f;

			// when p(z_n > 0)
			prior_i = q_ / ( float )LW1;								// p(z_n = i), i > 0
			for( i = 0; i < L; i++ ){									// i runs over all nucleotides in sequence
				y = posSeqs_[n].extractKmer( i, ( i < K_ ) ? i : K_ );	// extract (k+1)-mer y from positions (i-k,...,i)
				for( j = 0; j < ( ( W_ < i ) ? W_ : i ); j++ )			// j runs over all motif positions
					if( y != -1 )										// skip 'N' and other unknown alphabets
						r_[n][L-i+j] += s_[y][j];

					else
						break;
			}
			for( i = 1; i <= L; i++ ){
				r_[n][i] = expf( r_[n][i] ) * prior_i;
				normFactor += r_[n][i];
			}

			// when p(z_n = 0)
			r_[n][0] = prior_0;
			normFactor += r_[n][0];

			for( i = 0; i <= L; i++ )									// responsibility normalization
				r_[n][i] /= normFactor;

			llikelihood_ += logf( normFactor );
		}
	}

}

void EM::MStep(){

	int L, LW1, k, y, y2, j, n, i;

	// reset the fractional counts n
	for( k = 0; k < K_+1; k++ )
		for( y = 0; y < Global::powA[k+1]; y++ )
			for( j = 0; j < W_; j++ )
				n_[k][y][j] = 0.0f;

	// compute fractional occurrence counts for the highest order K
	if( Global::setSlow ){
		// slow code:
		for( n = 0; n < posSetN_; n++ ){								// n runs over all sequences
			L = posSeqs_[n].getL();
			LW1 = L - W_ + 1;
			for( i = 1; i <= LW1; i++ ){								// i runs over all nucleotides in sequence
				for( j = 0; j < W_; j++ ){								// j runs over all motif positions
					y = posSeqs_[n].extractKmer( i+j-1, ( j < K_ ) ? j : K_ );	// extract (k+1)-mer y
					if( y != -1 )										// skip 'N' and other unknown alphabets
						n_[K_][y][j] += r_[n][i];
					else
						break;
				}
			}
		}
	} else {
		// fast code:
		for( n = 0; n < posSetN_; n++ ){								// n runs over all sequences
			L = posSeqs_[n].getL();
			for( i = 0; i < L; i++ ){
				y = posSeqs_[n].extractKmer( i, ( i < K_ ) ? i : K_ );
				for( j = 0; j < ( ( W_ < i ) ? W_ : i ); j++ )
					if( y != -1 )										// skip 'N' and other unknown alphabets
						n_[K_][y][j] += r_[n][L-i+j];
					else
						break;
			}
		}
	}

	// compute fractional occurrence counts for lower orders
	for( k = K_; k > 0; k-- ){											// k runs over all orders
		for( y = 0; y < Global::powA[k+1]; y++ ){
			y2 = y % Global::powA[k];									// cut off the first nucleotide in (k+1)-mer
			for( j = 0; j < W_; j++ )
				n_[k-1][y2][j] += n_[k][y][j];
		}
	}

	// compute the Q function
//	Qfunc_ = calculateQfunc();

	// update model parameters v[k][y][j]
	motif_->updateV( n_, alpha_ );
}

void EM::optimizeAlphas(){
	// optimize alphas
	// motif.updateV()
}

void EM::optimizeQ(){
	// optimize Q function
	// motif.updateV()
}

float EM::calculateQfunc(){
	int n, L, LW1, i, j, y;
	float sumS = 0.0f, Qfunc = 0.0f;
	float prior_i, prior_0 = 1 - q_;

	for( n = 0; n < posSetN_; n++ ){
		L = posSeqs_[n].getL();
		LW1 = L - W_ + 1;
		prior_i = q_ / ( float )LW1;
		for( i = 0; i <= LW1; i++ ){
			for( j = 0; j < W_; j++ ){
				y =  posSeqs_[n].extractKmer( i, ( i < K_ ) ? i: K_ );
				sumS += s_[y][j];
			}
			Qfunc += r_[n][i] * sumS;
		}
		Qfunc += ( r_[n][0] * logf( prior_0 ) + ( 1 - r_[n][0] ) * logf( prior_i ) );
	}
	return Qfunc;
}

void EM::print(){
	for( int j = 0; j < W_; j++ ){
		for( int k = 0; k < K_+1; k++ ){
			for( int y = 0; y < Global::powA[k+1]; y++ )
				std::cout << std::scientific << v_motif_[k][y][j] << '\t';
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
}

void EM::write(){

	/**
	 * save EM parameters in five flat files:
	 * (1) posSequenceBasename.conds: 		conditional probabilities
	 * (2) posSequenceBasename.probs: 		probabilities
	 * (3) posSequenceBasename.EMcounts:	refined counts of (k+1)-mers
	 * (4) posSequenceBasename.EMposterior: responsibilities, posterior distributions
	 * (5) posSequenceBasename.EMalpha:		hyper-parameter alphas
	 * (6) posSequenceBasename.EMlogOdds:	log odds scores
	 */

	int k, y, j, n, i, L;

	std::string opath = std::string( Global::outputDirectory ) + '/'
			+ std::string( Global::posSequenceBasename );

	// output conditional probabilities v[k][y][j] and (k+1)-mer counts n[k][y][j]
	std::string opath_v = opath + ".conds";
	std::string opath_p = opath + ".probs";
	std::string opath_n = opath + ".EMcounts";
	std::ofstream ofile_v( opath_v.c_str() );
	std::ofstream ofile_p( opath_p.c_str() );
	std::ofstream ofile_n( opath_n.c_str() );
	for( j = 0; j < W_; j++ ){
		for( k = 0; k < K_+1; k++ ){
			for( y = 0; y < Global::powA[k+1]; y++ ){
				ofile_v << std::scientific << std::setprecision(8) << v_motif_[k][y][j] << ' ';
				ofile_p << std::scientific << std::setprecision(8) << p_[k][y][j] << ' ';
				ofile_n << std::scientific << n_[k][y][j] << ' ';
			}
			ofile_v << std::endl;
			ofile_p << std::endl;
			ofile_n << std::endl;
		}
		ofile_v << std::endl;
		ofile_p << std::endl;
		ofile_n << std::endl;
	}

	// output responsibilities r[n][i]
	std::string opath_r = opath + ".EMposterior";
	std::ofstream ofile_r( opath_r.c_str() );
	for( n = 0; n < posSetN_; n++ ){
		L = posSeqs_[n].getL();
		for( i = 0; i <= L; i++ )
			ofile_r << std::scientific << r_[n][i] << ' ';
		ofile_r << std::endl;
	}

	// output parameter alphas alpha[k][j]
	std::string opath_alpha = opath + ".EMalpha";
	std::ofstream ofile_alpha( opath_alpha.c_str() );
	for( k = 0; k < K_+1; k++ ){
		ofile_alpha << "k = " << k << std::endl;
		for( j = 0; j < W_; j++ )
			ofile_alpha << std::scientific << alpha_[k][j] << ' ';
		ofile_alpha << std::endl;
	}

	// output log odds scores s[y][j]
	std::string opath_s = opath + ".EMlogOdds";
	std::ofstream ofile_s( opath_s.c_str() );
	for( y = 0; y < Global::powA[K_+1]; y++ ){
		for( j = 0; j < W_; j++ )
			ofile_s << std::fixed << std::setprecision(6) << s_[y][j] << ' ';
		ofile_s << std::endl;
	}
}
