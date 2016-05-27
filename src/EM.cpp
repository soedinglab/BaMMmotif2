/*
 * EM.cpp
 *
 *  Created on: Apr 18, 2016
 *      Author: administrator
 */

#include "EM.h"

EM::EM( Motif* motif, BackgroundModel bg, std::vector<int> folds = std::vector<int> () ){
	motif_ = motif;
	bg_ = bg;
	if( folds.size() == 0 )		// for full set
		folds_ = std::vector<int> ( Global::cvFold );
	else						// for cross-validation
		folds_ = folds;

	powA_ = new int[Global::modelOrder+2];
	for( int k = 0; k < Global::modelOrder + 2; k++ )
		powA_[k] = Global::ipow( Alphabet::getSize(), k );

	Y_ = 0;
	for( int k = 0; k <= Global::modelOrder; k++ )
		Y_ += powA_[k+1]; 			// Y runs over all k-mers

	// allocate memory for s_[y][j] and initialize it
	s_ = ( float** )calloc( Y_, sizeof( float* ) );
	for( int y = 0; y < Y_; y++ ){
		s_[y] = ( float* )calloc( motif_->getW(), sizeof( float ) );
		for( int j = 0; j < motif_->getW(); j++ )
			s_[y][j] = 0.0f;
	}

    // allocate memory for r_[n][i] and initialize it
	r_ = ( float** )calloc( Global::posSequenceSet->getN(), sizeof( float* ) );
    for( unsigned int n = 0; n < Global::posSequenceSet->getN(); n++ ){
    	r_[n] = ( float* )calloc( Global::posSequenceSet->getSequences()[n].getL(), sizeof( float ) );
    	for( int i = 0; i < Global::posSequenceSet->getSequences()[n].getL(); i++ )
    		r_[n][i] = 0.0f;
    }

    // allocate memory for n_[k][y][j] and initialize it
	n_ = ( float*** )calloc( ( Global::modelOrder+1 ), sizeof( float** ) );
	for( int k = 0; k <= Global::modelOrder; k++ ){
		n_[k] = ( float** )calloc( powA_[k+1], sizeof( float* ) );
		for( int y = 0; y < powA_[k+1]; y++ ){
			n_[k][y] = ( float* )calloc( motif_->getW(), sizeof( float ) );
			for( int j = 0; j <  motif_->getW(); j++ )
				n_[k][y][j] = 0.0f;
		}
	}

	// allocate memory for alpha_[k][j] and initialize it
	alpha_ = ( float** )calloc( ( Global::modelOrder+1 ), sizeof( float* ) );
	for( int k = 0; k <= Global::modelOrder; k++ ){
		alpha_[k] = ( float* )calloc( motif_->getW(), sizeof( float ) );
		for( int j = 0; j < motif_->getW(); j++ )
			alpha_[k][j] = 0.0f;
	}

    // set initial values for alpha (alpha_k = beta x gamma^(order-1) with beta=20 and gamma=3) and q (0.9)
    q_ = 0.9f;
    likelihood_ = 0;		// default value?
    EMIterations_ = 0;		// default value?
}

EM::~EM(){				// Not the correct way to do that!!
	free( s_ );
	free( r_ );
	free( n_ );
	free( alpha_ );
}

int EM::learnMotif(){

	printf( " _____________________________\n"
			"|                             |\n"
			"| Compute log odds scores...  |\n"
			"|_____________________________|\n\n" );
	score();	// calculate s_[y][j]

    // iterate over
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
	if( !Global::noAlphaOptimization )
		optimizeAlphas();

    // * optional: optimizeQ()
	if( !Global::noQOptimization )
		optimizeQ();

    // * check likelihood for convergence


    // * remark: iterate over sequences using Global::posFoldIndices and folds_
	for( unsigned  int f=0; f < folds_.size() ; f++ ){
		int fold = folds_[f];
		for( unsigned int n=0; n < Global::posFoldIndices[fold].size(); n++ ){
			Sequence seq = Global::posSequenceSet->getSequences()[Global::posFoldIndices[fold][n]];
		}
	}

	if( Global::verbose ) print();

    return 0;
}

void EM::EStep(){

	// slow code to calculate responsibilities r_[n][i] at position i in sequence n

	// fast code to calculate responsibilities r_[n][i] at position i in sequence n
	int L, W, y;
	float p;
	for( int k = 0; k < Global::modelOrder; k++ ){
		for( unsigned int n = 0; n < Global::posSequenceSet->getN(); n++ ){	// n runs over all sequences
			Sequence seq = Global::posSequenceSet->getSequences()[n];
			L = seq.getL();
			W = motif_->getW();
			p = q_ / ( L - W + 1 );											// p(z_n = i)
			for( int i = k; i < L; i++ ){									// i runs over all nucleotides in sequence
				y = seq.extractKmer( i, k );								// extract (k+1)-mer y
				for( int j = 0; j < i; j++ )								// j runs over all motif positions
					r_[n][L-W+1-i+j] += s_[y][j];
				r_[n][L-i] = exp( r_[n][L-i] ) * p;
			}
			float sum = 0.0f;												// normalize r_[n][i]
			for( int i = 0; i < L-W; i++ ) sum += r_[n][i];
			for( int i = 0; i < L-W; i++ ) r_[n][i] /= sum;
		}
	}

}

void EM::MStep(){

	// resetN()
	int L, W, y;
	for( int k = 0; k < Global::modelOrder; k++ ){
		for( unsigned int n = 0; n < Global::posSequenceSet->getN(); n++ ){	// n runs over all sequences
			Sequence seq = Global::posSequenceSet->getSequences()[n];
			L = seq.getL();
			W = motif_->getW();
			for( int i = k; i < L; i++ ){									// i runs over all nucleotides in sequence
				y = seq.extractKmer( i, k );								// extract (k+1)-mer y
				for( int j = 0; j < i; j++ )								// j runs over all motif positions
					n_[k][y][j] += r_[n][L-W+1-i+j];
			}
		}
	}

	// M step


}

void EM::score(){

	// implementation of s[y][j]
	for( int k = 0; k <= Global::modelOrder; k++ ){
		std::cout << std::endl << "when k = " << k << std::endl;
		for( int y = 0; y < powA_[k+1]; y++ ){
			int y_bg = y - 1;
			for( int i = 0; i < k + 1; i++ )
				y_bg += powA_[i];
			for( int j = k; j < motif_->getW(); j++ ){
				s_[y][j] = log( motif_->getV()[k][y][j] / bg_.getV()[y_bg] );
				std::cout << s_[y][j] << '\t';
			}
			std::cout << std::endl;
		}
	}
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

