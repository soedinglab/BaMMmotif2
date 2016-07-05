/*
 * EM.cpp
 *
 *  Created on: Apr 18, 2016
 *      Author: administrator
 */

#include "EM.h"

EM::EM( Motif* motif, BackgroundModel bg, std::vector<int> folds = std::vector<int> () ){
	motif_ = new Motif( *motif );		// deep copy

	if( folds.size() != 0 )
		folds_ = folds;					// for cross-validation

	v_bg_ = ( float** )calloc( Global::bgModelOrder+1, sizeof( float* ) );
	for( int k = 0; k < Global::bgModelOrder+1; k++ ){			// deep copy
		v_bg_[k] = ( float* )calloc( Global::powA[k+1], sizeof( float ) );
		for( int y = 0; y < Global::powA[k+1]; y++ )
			v_bg_[k][y] = bg.getV()[k][y];
	}

	s_ = ( float** )calloc( Global::powA[Global::modelOrder+1], sizeof( float* ) );
	for( int y = 0; y < Global::powA[Global::modelOrder+1]; y++ )
		s_[y] = ( float* )calloc( motif_->getW(), sizeof( float ) );

    // allocate memory for r_[n][i] and initialize it
	r_ = ( float** )calloc( Global::posSequenceSet->getN(), sizeof( float* ) );
    for( unsigned int n = 0; n < Global::posSequenceSet->getN(); n++ )
    	r_[n] = ( float* )calloc( Global::posSequenceSet->getSequences()[n].getL() + 1, sizeof( float ) );

    // allocate memory for n_[k][y][j] and initialize it
	n_ = ( float*** )calloc( Global::modelOrder+1, sizeof( float** ) );
	for( int k = 0; k <= Global::modelOrder; k++ ){
		n_[k] = ( float** )calloc( Global::powA[k+1], sizeof( float* ) );
		for( int y = 0; y < Global::powA[k+1]; y++ )
			n_[k][y] = ( float* )calloc( motif_->getW(), sizeof( float ) );
	}

	// allocate memory for alpha_[k][j] and initialize it
	alpha_ = ( float** )malloc( ( Global::modelOrder+1 ) * sizeof( float* ) );
	for( int k = 0; k <= Global::modelOrder; k++ ){
		alpha_[k] = ( float* )malloc( motif_->getW() * sizeof( float ) );
		for( int j = 0; j < motif_->getW(); j++ )				// initialize alpha_ with global alpha
			alpha_[k][j] = Global::modelAlpha[k];
	}
}

EM::~EM(){
	delete[] motif_;

	for( int k = 0; k <= Global::bgModelOrder; k++ )
		free( v_bg_[k] );
	free( v_bg_ );

	for( int y = 0; y < Global::powA[Global::modelOrder+1]; y++ )
		free( s_[y] );
	free( s_ );

	for( unsigned int n = 0; n < Global::posSequenceSet->getN(); n++ )
		free( r_[n] );
	free( r_ );

	for( int k = 0; k <= Global::modelOrder; k++ ){
		for( int y = 0; y < Global::powA[k+1]; y++ )
			free( n_[k][y] );
		free( n_[k] );
	}
	free( n_ );

	for( int k = 0; k <= Global::modelOrder; k++ )
		free( alpha_[k] );
	free( alpha_ );
}

int EM::learnMotif(){

	bool iterate = true;										// flag for iterating before convergence
	float difference;											// model parameter difference for each EM iteration
	int K = Global::modelOrder;
	int W = motif_->getW();
	float** v_prior;											// hold the parameters before EM
	float** v_post;												// hold the parameters after EM

	// allocate for prior and posterior parameters v with the highest order
	v_prior = ( float** )calloc( Global::powA[K], sizeof( float* ) );
	v_post = ( float** )calloc( Global::powA[K], sizeof( float* ) );
	for( int y = 0; y < Global::powA[K]; y++ ){
		v_prior[y] = ( float* )calloc( W, sizeof( float ) );
		v_post[y] = ( float* )calloc( W, sizeof( float ) );
	}

	while( iterate && ( EMIterations_ < Global::maxEMIterations ) ){
		// iterate over
		EMIterations_++;										// count times of iteration

		// before EM, get the prior parameter variables
		for( int y = 0; y < Global::powA[K]; y++ )
			for( int j = 0; j < W; j++ )
				v_prior[y][j] = motif_->getV()[K][y][j];

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
		for( int y = 0; y < Global::powA[K]; y++ )
			for( int j = 0; j < W; j++ )
				v_post[y][j] = motif_->getV()[K][y][j];


		// * check likelihood for convergence
		printf( "\n*********** Check convergence *************\n" );

		difference = 0.0f;										// reset difference to 0
		for( int y = 0; y < Global::powA[K]; y++ )
			for( int j = 0; j < W; j++ )
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

	if( Global::verbose ) print();
    return 0;
}

void EM::EStep(){

	int L;																// length of sequences
	int W = motif_->getW();
	int K = Global::modelOrder;
	float prior_i;														// positional preference prior state for p(z_n = i), i > 0
	float prior_0 = 1 - q_;												// positional preference prior state for p(z_n = 0), no motif present

	// compute log odd scores s[y][j]
	for( int y = 0; y < Global::powA[K+1]; y++ ){
		for( int j = 0; j < W; j++ )
			if( K <= Global::bgModelOrder ){
				s_[y][j] = log( motif_->getV()[K][y][j] / v_bg_[K][y] );
			} else {
				int y3 = y % Global::powA[3];							// 3 rightmost nucleotides in (k+1)-mer y
				s_[y][j] = log( motif_->getV()[K][y][j] / v_bg_[Global::bgModelOrder][y3] );
			}
	}

	// calculate responsibilities r_[n][i] at position i in sequence n
//	// slow code:
//	printf( "\n*********** Slow code: E-Step ***********\n" );
//	for( unsigned int n = 0; n < Global::posSequenceSet->getN(); n++ ){	// n runs over all sequences
//		Sequence seq = Global::posSequenceSet->getSequences()[n];
//		L = seq.getL();
//
//		// when p(z_n > 0)
//		prior_i = q_ / ( L - W + 1 );									// p(z_n = i), i > 0
//		int y;
//		float sum = 0.0f;												// normalize r_[n][i]
//		for( int i = 0; i < L-W+1; i++ ){
//			for( int j = 0; j < W; j++ ){
//				y = seq.extractKmer( i+j, K );
//				if( y != -1 )											// skip 'N' and other unknown alphabets
//					r_[n][i] += s_[y][j];
//				else
//					break;
//			}
//			r_[n][i] = exp( r_[n][i] ) * prior_i;
//			sum += r_[n][i];
//		}
//		// when p(z_n = 0)
//		r_[n][L] = prior_0;
//		sum += r_[n][L];
//
//		for( int i = 0; i < L*-W+1; i++ ){
//			r_[n][i] /= sum;
//		}
//	}

	// fast code:
	printf( "\n*********** Fast code: E-Step ***********\n" );
	for( unsigned int n = 0; n < Global::posSequenceSet->getN(); n++ ){	// n runs over all sequences
		Sequence seq = Global::posSequenceSet->getSequences()[n];
		L = seq.getL();

		// when p(z_n > 0)
		prior_i = q_ / ( L - W + 1 );
		int y;
		for( int i = 0; i < L; i++ ){									// i runs over all nucleotides in sequence
			y = seq.extractKmer( i, K );								// extract (k+1)-mer y from positions (i-k,...,i)
			if( y != -1 )												// skip 'N' and other unknown alphabets
				for( int j = 0; j < W; j++ )							// j runs over all motif positions
					r_[n][L-W-i+j] += s_[y][j];
			else
 				break;
		}
		float sum = 0.0f;
		for( int i = 0; i < L; i++ ){
			r_[n][i] = exp( r_[n][i] ) * prior_i;
			sum += r_[n][i];
		}
		// when p(z_n = 0)
		r_[n][L] = prior_0;
		sum += r_[n][L];

		for( int i = 0; i <= L; i++ )									// normalization
			r_[n][i] /= sum;
	}

}

void EM::MStep(){
	int L;
	int W = motif_->getW();
	int K = Global::modelOrder;

	// reset the fractional counts n
	for( int k = 0; k <= K; k++ )
		for( int y = 0; y < Global::powA[k+1]; y++ )
			for( int j = 0; j < W; j++ )
				n_[k][y][j] = 0.0f;

	// compute fractional occurrence counts for the highest order K
//	// slow code:
//	printf( "\n*********** Slow code: M-Step ***********\n" );
//	for( unsigned int n = 0; n < Global::posSequenceSet->getN(); n++ ){	// n runs over all sequences
//		Sequence seq = Global::posSequenceSet->getSequences()[n];
//		L = seq.getL();
//		int y;
//		for( int i = 0; i < L-W+1; i++ ){								// i runs over all nucleotides in sequence
//			for( int j = 0; j < W; j++ ){								// j runs over all motif positions
//				y = seq.extractKmer( i+j, K );							// extract (k+1)-mer y
//				if( y != -1 )											// skip 'N' and other unknown alphabets
//					n_[K][y][j] += r_[n][i];
//				else
//					break;
//			}
//		}
//	}

	// fast code:
	printf( "\n*********** Fast code: M-Step ***********\n" );
	for( unsigned int n = 0; n < Global::posSequenceSet->getN(); n++ ){	// n runs over all sequences
		Sequence seq = Global::posSequenceSet->getSequences()[n];
		L = seq.getL();
		int y;
		for( int i = 0; i < L; i++ ){
			y = seq.extractKmer( i, K );
			for( int j = 0; j < W; j++ ){
				if( y != -1 )											// skip 'N' and other unknown alphabets
					n_[K][y][j] += r_[n][L-W-i+j];
				else
					break;
			}
		}
	}

	// compute fractional occurrence counts for lower orders
	for( int k = K; k > 0; k-- ){										// k runs over all orders
		for( int y = 0; y < Global::powA[k+1]; y++ ){
			int y2 = y % Global::powA[k];								// cut off the first nucleotide in (k+1)-mer
			for( int j = 0; j < W; j++ )
				n_[k-1][y2][j] += n_[k][y][j];
			int yk = y / Global::powA[1];								// cut off the last nucleotide in (k+1)-mer
			n_[k-1][yk][-1] += n_[k][y][0];
		}
	}

	// only for testing:
	for(int k=0; k <= K; k++)
		for(int y = 0; y < Global::powA[k+1]; y++)
			for(int j = 0; j < W; j++){
				if( n_[k][y][j] < 0 )
				std::cout << "n_["<<k<<"]["<<y<<"]["<<j<<"] = " << n_[k][y][j] << '\t';
			}

	// update model parameters v[k][y][j]
	motif_->updateV( n_, alpha_ );

}

void EM::optimizeAlphas(){
	// optimize alphas
	// motif.updateV()
}

void EM::optimizeQ(){

}

void EM::print(){

}

void EM::write(){

}
