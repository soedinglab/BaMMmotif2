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

	powA_ = new int[Global::modelOrder+2];
	for( int k = 0; k < Global::modelOrder + 2; k++ )
		powA_[k] = Global::ipow( Alphabet::getSize(), k );

	Y_ = 0;
	for( int k = 0; k <= Global::modelOrder; k++ )
		Y_ += powA_[k+1]; 				// Y runs over all k-mers

	v_bg_ = ( float* )calloc( Y_, sizeof( float ) );
	for( int y = 0; y< Y_; y++ )		// deep copy
		v_bg_[y] = bg.getV()[y];

    // allocate memory for r_[n][i] and initialize it
	r_ = ( float** )calloc( Global::posSequenceSet->getN(), sizeof( float* ) );
    for( unsigned int n = 0; n < Global::posSequenceSet->getN(); n++ ){
    	r_[n] = ( float* )calloc( Global::posSequenceSet->getSequences()[n].getL(), sizeof( float ) );
    }

    // allocate memory for n_[k][y][j] and initialize it
	n_ = ( float*** )calloc( ( Global::modelOrder+1 ), sizeof( float** ) );
	for( int k = 0; k <= Global::modelOrder; k++ ){
		n_[k] = ( float** )calloc( powA_[k+1], sizeof( float* ) );
		for( int y = 0; y < powA_[k+1]; y++ ){
			n_[k][y] = ( float* )calloc( motif_->getW(), sizeof( float ) );
		}
	}

	// allocate memory for alpha_[k][j] and initialize it
	alpha_ = ( float** )malloc( ( Global::modelOrder+1 ) * sizeof( float* ) );
	for( int k = 0; k <= Global::modelOrder; k++ ){
		alpha_[k] = ( float* )malloc( motif_->getW() * sizeof( float ) );
		for( int j = 0; j < motif_->getW(); j++ )			// initialize alpha_ with global setting
			alpha_[k][j] = Global::modelAlpha[k];
	}

    // set initial values for alpha (alpha_k = beta x gamma^(order-1) with beta=20 and gamma=3) and q (0.9)
    q_ = 0.9f;
    likelihood_ = 0;		// default value?
//    EMIterations_ = 0;		// default value?
}

EM::~EM(){
	delete[] motif_;

	free( v_bg_ );

	for( unsigned int n = 0; n < Global::posSequenceSet->getN(); n++ )
		free( r_[n] );
	free( r_ );

	for( int k = 0; k <= Global::modelOrder; k++ ){
		for( int y = 0; y < powA_[k+1]; y++ )
			free( n_[k][y] );
		free( n_[k] );
	}
	free( n_ );

	for( int k = 0; k <= Global::modelOrder; k++ )
		free( alpha_[k] );
	free( alpha_ );

	delete[] powA_;
}

int EM::learnMotif(){

	bool iterate = true;			// flag for iterating before convergence
	unsigned int EMIterations = 0;	// count times of iteration
	float difference;				// model parameter difference between EM iterations

	float sumV_prior, sumV_post;

	while( iterate && ( Global::maxEMIterations > EMIterations ) ){

		sumV_prior = sumV( motif_ );		// hold the sum of v before EM

		// iterate over
		EMIterations++;				// count times of iteration
		printf( " _____________________________\n"
				"|                             |\n"
				"| E-step: calculate posterior |\n"
				"|_____________________________|\n\n" );
		EStep();

		printf( " ___________________________\n"
				"|                           |\n"
				"| M-step: update parameters |\n"
				"|___________________________|\n\n" );
		MStep();


		// * optional: optimizeAlphas()
		if( !Global::noAlphaOptimization ){
			printf( " ____________________\n"
					"|                    |\n"
					"| Optimize alphas... |\n"
					"|____________________|\n\n" );
			optimizeAlphas();
		}

		// * optional: optimizeQ()
		if( !Global::noQOptimization ){
			printf( " ________________________\n"
					"|                        |\n"
					"| Optimize Q function... |\n"
					"|________________________|\n\n" );
			optimizeQ();
		}

		sumV_post = sumV( motif_ );			// hold the sum of v after EM

		difference = sumV_post - sumV_prior;
		std::cout << "difference = " << difference << std::endl;
		// * check likelihood for convergence
		printf( " ______________________\n"
				"|                      |\n"
				"| Check convergence... |\n"
				"|______________________|\n\n" );
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

	// slow code to calculate responsibilities r_[n][i] at position i in sequence n
	int L, y;
	int W = motif_->getW();
	float p;
	for( unsigned int n = 0; n < Global::posSequenceSet->getN(); n++ ){	// n runs over all sequences
		Sequence seq = Global::posSequenceSet->getSequences()[n];
		L = seq.getL();
		p = q_ / ( L - W + 1 );											// p(z_n = i)
		for( int i = Global::modelOrder; i < L; i++ ){
			y = seq.extractKmer( i, Global::modelOrder );
			if( y != -1 ){												// skip 'N' and other unknown alphabets
				for( int j = 0; j < W; j++ )
					r_[n][i-j] += motif_->getS()[y][j];
				r_[n][i] = exp( r_[n][i] ) * p;
			}
			else
				continue;
		}

		float sum = 0.0f;												// normalize r_[n][i]
		for( int i = 0; i < L-W; i++ ) sum += r_[n][i];
		for( int i = 0; i < L-W; i++ ) r_[n][i] /= sum;
//		for( int i = 0; i < L-W; i++ ) std::cout << r_[n][i] << '\t';
//		std::cout << std::endl << "when n = " << n << ", i = " << 3 << ", r[n][3] = " << r_[n][3] << std::endl;
	}

//	// fast code to calculate responsibilities r_[n][i] at position i in sequence n
//	int L, y;
//	int W = motif_->getW();
//	float p;
//	for( unsigned int n = 0; n < Global::posSequenceSet->getN(); n++ ){	// n runs over all sequences
//		Sequence seq = Global::posSequenceSet->getSequences()[n];
//		L = seq.getL();
//		p = q_ / ( L - W + 1 );											// p(z_n = i)
//		for( int i = Global::modelOrder; i < L; i++ ){					// i runs over all nucleotides in sequence
//			y = seq.extractKmer( i, Global::modelOrder );				// extract (k+1)-mer y
//			if( y != -1 ){												// skip 'N' and other unknown alphabets
//				for( int j = 0; j < i; j++ )							// j runs over all motif positions
//					r_[n][L-W+1-i+j] += motif_->getS()[y][j];
//				r_[n][L-i] = exp( r_[n][L-i] ) * p;
//			}
//			else
//				continue;
//		}
//		float sum = 0.0f;												// normalize r_[n][i]
//		for( int i = 0; i < L-W; i++ ) sum += r_[n][i];
//		for( int i = 0; i < L-W; i++ ) r_[n][i] /= sum;
//		for( int i = 0; i < L-W; i++ ) std::cout << r_[n][i] << '\t';
//		std::cout << std::endl << "when n = " << n << ", i = " << 3 << ", r[n][3] = " << r_[n][3] << std::endl;
//	}

}

void EM::MStep(){

	// reset n, the fractional counts
	printf( " _________________\n"
			"|                 |\n"
			"| resetting n_occ |\n"
			"|_________________|\n\n" );
	int L, y;
	int W = motif_->getW();
	int K = Global::modelOrder;
	// compute fractional occurrence counts for the highest order K
	for( unsigned int n = 0; n < Global::posSequenceSet->getN(); n++ ){	// n runs over all sequences
		Sequence seq = Global::posSequenceSet->getSequences()[n];
		L = seq.getL();
		for( int i = K; i < L; i++ ){									// i runs over all nucleotides in sequence
			y = seq.extractKmer( i, K );								// extract (k+1)-mer y
			if( y != -1 ){												// skip 'N' and other unknown alphabets
				for( int j = 0; j < i; j++ )							// j runs over all motif positions
					n_[K][y][j] += r_[n][L-W+1-i+j];
			} else
				continue;
		}
	}

	// compute fractional occurrence counts for lower orders
	for( int k = K; k > 0; k-- )										// k runs over all orders
		for( int y = 0; y < powA_[k+1]; y++ ){
			int y2 = y / powA_[1];										// cut off the first nucleotide in (k+1)-mer
			for( int j = k-1; j < W; j++ )
				n_[k-1][y2][j] += n_[k][y][j];
			int yk = y % powA_[k];										// cut off the last nucleotide in (k+1)-mer
			n_[k-1][yk][-1] += n_[k][y][0];
		}

	// update model parameters v[k][y][j]
	printf( " _______________\n"
			"|               |\n"
			"| updating v... |\n"
			"|_______________|\n\n" );

	// No idea why here alpha_ cannot be passed...
	float** alpha;
	alpha = ( float** )malloc( ( Global::modelOrder+1 ) * sizeof( float* ) );
	for( int k = 0; k <= Global::modelOrder; k++ ){
		alpha[k] = ( float* )malloc( motif_->getW() * sizeof( float ) );
		for( int j = 0; j < motif_->getW(); j++ )			// initialize alpha_ with global setting
			alpha[k][j] = Global::modelAlpha[k];
	}

	motif_->updateV( n_, alpha );
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

