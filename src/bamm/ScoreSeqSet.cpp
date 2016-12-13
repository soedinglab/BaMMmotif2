/*
 * ScoreSeqSet.cpp
 *
 *  Created on: Dec 13, 2016
 *      Author: wanwan
 */

#include "ScoreSeqSet.h"

ScoreSeqSet::ScoreSeqSet( Motif* motif, BackgroundModel* bg, std::vector<Sequence*> seqSet ){

	for( int k = 0; k < std::max( Global::modelOrder+2 , Global::bgModelOrder+2 ); k++ ){
		Y_.push_back( ipow( Alphabet::getSize(), k ) );
	}

	motif_ = motif;
	bg_ = bg;
	seqSet_ = seqSet;

	// allocate memory for scores_[n][i]
	scores_ = ( float** )calloc( seqSet_.size(), sizeof( float* ) );
	for( size_t n = 0; n < seqSet_.size(); n++ ){
		int LW1 = seqSet_[n]->getL() - motif_->getW() + 1;
		scores_[n] = ( float* )calloc( LW1, sizeof( float ) );
	}
}

ScoreSeqSet::~ScoreSeqSet(){

	for( size_t n = 0; n < seqSet_.size(); n++ ){
		free( scores_[n] );
	}
	free( scores_ );

}

// store the log odds scores at all positions of each sequence
float** ScoreSeqSet::score(){

	int K = Global::modelOrder;
	int K_bg = Global::bgModelOrder;
	int W = motif_->getW();

	for( size_t n = 0; n < seqSet_.size(); n++ ){

		int LW1 = seqSet_[n]->getL() - W + 1;

		for( int i = 0; i < LW1; i++ ){

			for( int j = 0; j < W; j++ ){

				int y = seqSet_[n]->extractKmer( i+j, std::min(i+j, K ) );

				int y_bg = y % Y_[K_bg+1];

				if( y >= 0 ){

					scores_[n][i] += ( logf( motif_->getV()[K][y][j] ) - logf( bg_->getV()[std::min( K, K_bg )][y_bg] ) );

				}
			}
		}

	}
	return scores_;
}

void ScoreSeqSet::write(){

	/**
	 * save log odds scores in one flat file:
	 * posSequenceBasename.logOdds
	 */

	std::string opath = std::string( Global::outputDirectory ) + '/'
			+ Global::posSequenceBasename + ".logOdds";

	std::ofstream ofile( opath );

	for( size_t n = 0; n < seqSet_.size(); n++ ){
		int LW1 = seqSet_[n]->getL() - motif_->getW() + 1;
		for( int i = 0; i < LW1; i++ ){
			ofile << std::setprecision( 3 ) << scores_[n][i] << '\t';
		}
		ofile << std::endl;
	}

}
