/*
 * EM.cpp
 *
 *  Created on: Apr 18, 2016
 *      Author: administrator
 */

#include "EM.h"

EM::EM( Motif* motif, BackgroundModel* bg, std::vector<int> folds ){
	motif_ = motif;												// deep copy
	bg_ = bg;

	W_ = motif_->getW();
	if( folds.size() != 0 )
		folds_ = folds;											// for cross-validation

//	v_bg_ = ( float** )calloc( k_bg_+1, sizeof( float* ) );
//	for( int k = 0; k < k_bg_+1; k++ ){							// deep copy
//		v_bg_[k] = ( float* )calloc( Global::powA[k+1], sizeof( float ) );
//		for( int y = 0; y < Global::powA[k+1]; y++ )
//			v_bg_[k][y] = bg->getVbg()[k][y];
//	}

	s_ = ( float** )calloc( Global::powA[K_+1], sizeof( float* ) );
	for( int y = 0; y < Global::powA[K_+1]; y++ )
		s_[y] = ( float* )calloc( W_, sizeof( float ) );

    // allocate memory for r_[n][i] and initialize it
	r_ = ( float** )calloc( Global::posSequenceSet->getN(), sizeof( float* ) );
    for( unsigned int n = 0; n < Global::posSequenceSet->getN(); n++ ){
    	r_[n] = ( float* )calloc( Global::posSequenceSet->getSequences()[n].getL() + 1/* + W_*/, sizeof( float ) );
    	for( int i = 0; i < Global::posSequenceSet->getSequences()[n].getL() + 1; i++ )
    		r_[n][i] = 0.0f;
    }

    // allocate memory for n_[k][y][j] and initialize it
	n_ = ( float*** )calloc( K_+1, sizeof( float** ) );
	for( int k = 0; k < K_+1; k++ ){
		n_[k] = ( float** )calloc( Global::powA[k+1], sizeof( float* ) );
		for( int y = 0; y < Global::powA[k+1]; y++ )
			n_[k][y] = ( float* )calloc( W_, sizeof( float ) );
	}

	// allocate memory for alpha_[k][j] and initialize it
	alpha_ = ( float** )malloc( ( K_+1 ) * sizeof( float* ) );
	for( int k = 0; k < K_+1; k++ ){
		alpha_[k] = ( float* )malloc( W_ * sizeof( float ) );
		for( int j = 0; j < W_; j++ )							// initialize alpha_ with global alpha
			alpha_[k][j] = Global::modelAlpha[k];
	}
}

EM::~EM(){

//	for( int k = 0; k <= k_bg_; k++ )
//		free( v_bg_[k] );
//	free( v_bg_ );

	for( int y = 0; y < Global::powA[K_+1]; y++ )
		free( s_[y] );
	free( s_ );

	for( unsigned int n = 0; n < Global::posSequenceSet->getN(); n++ )
		free( r_[n] );
	free( r_ );

	for( int k = 0; k <= K_; k++ ){
		for( int y = 0; y < Global::powA[k+1]; y++ )
			free( n_[k][y] );
		free( n_[k] );
	}
	free( n_ );

	for( int k = 0; k <= K_; k++ )
		free( alpha_[k] );
	free( alpha_ );
	std::cout << "Destructor for EM class works fine. \n";
}

int EM::learnMotif(){

	bool iterate = true;										// flag for iterating before convergence
	float difference;											// model parameter difference for each EM iteration
	float** v_prior;											// hold the parameters before EM
	float** v_post;												// hold the parameters after EM

	// allocate for prior and posterior parameters v with the highest order
	v_prior = ( float** )calloc( Global::powA[K_+1], sizeof( float* ) );
	v_post = ( float** )calloc( Global::powA[K_+1], sizeof( float* ) );
	for( int y = 0; y < Global::powA[K_+1]; y++ ){
		v_prior[y] = ( float* )calloc( W_, sizeof( float ) );
		v_post[y] = ( float* )calloc( W_, sizeof( float ) );
	}

	while( iterate && ( EMIterations_ < Global::maxEMIterations ) ){
		// iterate over
		EMIterations_++;										// count times of iteration

		// before EM, get the prior parameter variables
		for( int y = 0; y < Global::powA[K_+1]; y++ )
			for( int j = 0; j < W_; j++ )
				v_prior[y][j] = motif_->getV()[K_][y][j];

		printf( " ______\n"
				"|      |\n"
				"|  EM  |\n"
				"|______|\n\n" );

		// E-step: calculate posterior
		EStep();

		// M-step: update parameters
		MStep();

		// * optional: optimizeAlphas()
		if( !Global::noAlphaOptimization )	optimizeAlphas();

		// * optional: optimizeQ()
		if( !Global::noQOptimization )		optimizeQ();


		// after EM, get the posterior parameter variables
		for( int y = 0; y < Global::powA[K_+1]; y++ )
			for( int j = 0; j < W_; j++ )
				v_post[y][j] = motif_->getV()[K_][y][j];

		// * check likelihood for convergence
		printf( "\n*********** Check convergence ***********\n" );

		difference = 0.0f;										// reset difference to 0
		for( int y = 0; y < Global::powA[K_+1]; y++ )
			for( int j = 0; j < W_; j++ )
				difference += fabs( v_post[y][j] - v_prior[y][j] );
		std::cout << "difference = " << difference <<  " at the " << EMIterations_ << "th iteration." << std::endl;
		if( difference < Global::epsilon )	iterate = false;
	}

    // * remark: iterate over sequences using Global::posFoldIndices and folds_
//	if( folds_.size() > 0 )
//		for( unsigned  int f = 0; f < folds_.size() ; f++ ){
//			int fold = folds_[f];
//			for( unsigned int n = 0; n < Global::posFoldIndices[fold].size(); n++ )
//				Sequence seq = Global::posSequenceSet->getSequences()[Global::posFoldIndices[fold][n]];
//		}

//	if( Global::verbose ) print();

	for( int y = 0; y < Global::powA[K_+1]; y++ ){
		free( v_prior[y] );
		free( v_post[y] );
	}
	free( v_prior );
	free( v_post );

    return 0;
}

void EM::EStep(){

	int L;																// length of sequences
	float prior_i;														// positional preference prior state for p(z_n = i), i > 0
	float prior_0 = 1 - q_;												// positional preference prior state for p(z_n = 0), no motif present
	float normFactor;													// normalize responsibilities r[n][i]

	// compute log odd scores s[y][j], log likelihoods
	for( int y = 0; y < Global::powA[K_+1]; y++ ){
		for( int j = 0; j < W_; j++ ){
			if( K_ <= k_bg_ ){
//				s_[y][j] = log( motif_->getV()[K_][y][j] / v_bg_[K_][y] );
				s_[y][j] = log( motif_->getV()[K_][y][j] / bg_->getVbg()[K_][y] );
			} else {
				int y3 = y % Global::powA[3];							// 3 rightmost nucleotides in (k+1)-mer y
//				s_[y][j] = log( motif_->getV()[K_][y][j] / v_bg_[k_bg_][y3] );
				s_[y][j] = log( motif_->getV()[K_][y][j] / bg_->getVbg()[k_bg_][y3] );
			}
		}
	}

	// calculate responsibilities r_[n][i] at position i in sequence n
	// slow code:
	printf( "\n*********** Slow code: E-Step ***********\n" );
	for( unsigned int n = 0; n < Global::posSequenceSet->getN(); n++ ){	// n runs over all sequences
		Sequence seq = Global::posSequenceSet->getSequences()[n];
		L = seq.getL();

		// when p(z_n > 0)
		int y;
		prior_i = q_ / ( L - W_ + 1 );									// p(z_n = i), i > 0
		normFactor = 0.0f;												// reset normalization factor to 0
		for( int i = 0; i < L-W_+1; i++ ){
			for( int j = 0; j < W_; j++ ){
				y = seq.extractKmer( i+j, ( j < K_ ) ? j : K_ );
				if( y != -1 )											// skip 'N' and other unknown alphabets
					r_[n][i] += s_[y][j];
				else
					break;
			}
			r_[n][i] = exp( r_[n][i] ) * prior_i;
			normFactor += r_[n][i];
		}
		// when p(z_n = 0)
		r_[n][L-W_+1] = prior_0;
		normFactor += r_[n][L];

		for( int i = 0; i <= L-W_+1; i++ )								// normalization
			r_[n][i] /= normFactor;
	}

	// fast code:
//	printf( "\n*********** Fast code: E-Step ***********\n" );
//	for( unsigned int n = 0; n < Global::posSequenceSet->getN(); n++ ){	// n runs over all sequences
//		Sequence seq = Global::posSequenceSet->getSequences()[n];
//		L = seq.getL();
//
//		// when p(z_n > 0)
//		prior_i = q_ / ( L - W_ + 1 );
//		int y;
//		for( int i = 0; i < L; i++ ){									// i runs over all nucleotides in sequence
//			y = seq.extractKmer( i, K_ );								// extract (k+1)-mer y from positions (i-k,...,i)
//			if( y != -1 )												// skip 'N' and other unknown alphabets
//				for( int j = K_; j < W_; j++ )							// j runs over all motif positions
//					r_[n][L-W_-i+j] += s_[y][j];
//			else
// 				break;
//		}
//		normFactor = 0.0f;
//		for( int i = 0; i < L; i++ ){
//			r_[n][i] = exp( r_[n][i] ) * prior_i;
//			normFactor += r_[n][i];
//		}
//		// when p(z_n = 0)
//		r_[n][L] = prior_0;
//		normFactor += r_[n][L];
//
//		for( int i = 0; i <= L; i++ )									// normalization
//			r_[n][i] /= normFactor;
//	}
}

void EM::MStep(){
	int L;

	// reset the fractional counts n
	for( int k = 0; k < K_+1; k++ )
		for( int y = 0; y < Global::powA[k+1]; y++ )
			for( int j = 0; j < W_; j++ )
				n_[k][y][j] = 0.0f;

	// compute fractional occurrence counts for the highest order K
	// slow code:
	printf( "\n*********** Slow code: M-Step ***********\n" );
	for( unsigned int n = 0; n < Global::posSequenceSet->getN(); n++ ){	// n runs over all sequences
		Sequence seq = Global::posSequenceSet->getSequences()[n];
		L = seq.getL();
		int y;
		for( int i = 0; i < L-W_+1; i++ ){								// i runs over all nucleotides in sequence
			for( int j = 0; j < W_; j++ ){								// j runs over all motif positions
				y = seq.extractKmer( i+j, ( j < K_ ) ? j : K_ );		// extract (k+1)-mer y
				if( y != -1 )											// skip 'N' and other unknown alphabets
					n_[K_][y][j] += r_[n][i];
				else
					break;
			}
		}
	}

	// fast code:
//	printf( "\n*********** Fast code: M-Step ***********\n" );
//	for( unsigned int n = 0; n < Global::posSequenceSet->getN(); n++ ){	// n runs over all sequences
//		Sequence seq = Global::posSequenceSet->getSequences()[n];
//		L = seq.getL();
//		int y;
//		for( int i = 0; i < L; i++ ){
//			y = seq.extractKmer( i, K_ );
//			for( int j = K_; j < W_; j++ ){
//				if( y != -1 )											// skip 'N' and other unknown alphabets
//					n_[K_][y][j] += r_[n][L-W_-i+j];
//				else
//					break;
//			}
//		}
//	}

	// compute fractional occurrence counts for lower orders
	for( int k = K_; k > 0; k-- ){										// k runs over all orders
		for( int y = 0; y < Global::powA[k+1]; y++ ){
			int y2 = y % Global::powA[k];								// cut off the first nucleotide in (k+1)-mer
//			int yk = y / Global::powA[1];								// cut off the last nucleotide in (k+1)-mer
			for( int j = 0; j < W_; j++ )
				n_[k-1][y2][j] += n_[k][y][j];
//			n_[k-1][yk][-1] += n_[k][y][0];
		}
	}

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

void EM::print(){
	for( int j = 0; j < W_; j++ ){
		for( int k = 0; k < K_+1; k++ ){
			for( int y = 0; y < Global::powA[k+1]; y++ )
				std::cout << std::scientific << motif_->getV()[k][y][j] << '\t';
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
}

void EM::write(){
	/*
	 * save different EM parameters in five flat files:
	 * (1) posSequenceBasename.EMconds: 	conditional probabilities
	 * (2) posSequenceBasename.EMcounts:	refined counts of (k+1)-mers
	 * (3) posSequenceBasename.EMresp: 		responsibilities
	 * (4) posSequenceBasename.EMalpha:		hyper-parameter alphas
	 * (5) posSequenceBasename.EMlogOdds:	log odds scores
	 */
	std::string opath = std::string( Global::outputDirectory ) + '/'
			+ std::string( Global::posSequenceBasename );

	// output conditional probabilities v[k][y][j] and n[k][y][j]
	std::string opath_v = opath + ".EMconds";
	std::string opath_n = opath + ".EMcounts";
	std::ofstream ofile_v( opath_v.c_str() );
	std::ofstream ofile_n( opath_n.c_str() );
	for( int j = 0; j < W_; j++ ){
		for( int k = 0; k < K_+1; k++ ){
			for( int y = 0; y < Global::powA[k+1]; y++ ){
				ofile_v << std::scientific << motif_->getV()[k][y][j] << '\t';
				ofile_n << std::scientific << n_[k][y][j] << '\t';
			}
			ofile_v << std::endl;
			ofile_n << std::endl;
		}
		ofile_v << std::endl;
		ofile_n << std::endl;
	}

	// output responsibilities r[n][i]
	std::string opath_r = opath + ".EMresp";
	std::ofstream ofile_r( opath_r.c_str() );
	for( unsigned int n = 0; n < Global::posSequenceSet->getN(); n++ ){
		Sequence seq = Global::posSequenceSet->getSequences()[n];
		int L = seq.getL();
		for( int i = 0; i < L-W_+1; i++ )
			ofile_r << std::scientific << r_[n][i] << '\t';
		ofile_r << std::endl;
	}

	// output parameter alphas alpha[k][j]
	std::string opath_alpha = opath + ".EMalpha";
	std::ofstream ofile_alpha( opath_alpha.c_str() );
	for( int k = 0; k < K_+1; k++ ){
		ofile_alpha << "k = " << k << std::endl;
		for( int j = 0; j < W_; j++ )
			ofile_alpha << std::scientific << alpha_[k][j] << '\t';
		ofile_alpha << std::endl;
	}

	// output log odds scores s[y][j]
	std::string opath_s = opath + ".EMlogOdds";
	std::ofstream ofile_s( opath_s.c_str() );
	for( int y = 0; y < Global::powA[K_+1]; y++ ){
		for( int j = 0; j < W_; j++ )
			ofile_s << std::fixed << std::setprecision(6) << s_[y][j] << '\t';
		ofile_s << std::endl;
	}
}
