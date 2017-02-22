/*
 * ScoreSeqSet.cpp
 *
 *  Created on: Dec 13, 2016
 *      Author: wanwan
 */

#include "ScoreSeqSet.h"

#include <float.h>		// -FLT_MAX

ScoreSeqSet::ScoreSeqSet( Motif* motif, BackgroundModel* bg, std::vector<Sequence*> seqSet ){

	for( int k = 0; k < std::max( Global::modelOrder+2 , Global::bgModelOrder+2 ); k++ ){
		Y_.push_back( ipow( Alphabet::getSize(), k ) );
	}

	motif_ = motif;
	bg_ = bg;
	seqSet_ = seqSet;

	// allocate memory for scores_[n][i]
	scores_.resize( 2 );

}

ScoreSeqSet::~ScoreSeqSet(){

}

// store the log odds scores at all positions of each sequence
void ScoreSeqSet::score(){

	int 	K = Global::modelOrder;
	int 	K_bg = Global::bgModelOrder;
	int 	W = motif_->getW();
	float 	maxScore;								// maximal logOddsScore over all positions for each sequence

	for( size_t n = 0; n < seqSet_.size(); n++ ){

		int LW1 = seqSet_[n]->getL() - W + 1;

		maxScore = -FLT_MAX;

		std::vector<float> logOdds( LW1 );

		for( int i = 0; i < LW1; i++ ){

			logOdds[i] = 0.0f;

			for( int j = 0; j < W; j++ ){

				int y = seqSet_[n]->extractKmer( i+j, std::min(i+j, K ) );

				int y_bg = y % Y_[K_bg+1];

				if( y >= 0 ){

					logOdds[i] += ( logf( motif_->getV()[K][y][j] ) - logf( bg_->getV()[std::min( K, K_bg )][y_bg] ) );

				}
			}

			// take all the log odds scores for MOPS model:
			scores_[0].push_back( logOdds[i] );

			// take the largest log odds score for ZOOPS model:
			if( logOdds[i] > maxScore ){
				maxScore = logOdds[i];
			}

		}

		scores_[1].push_back( maxScore );

	}
}

std::vector<std::vector<float>> ScoreSeqSet::getScores(){
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
