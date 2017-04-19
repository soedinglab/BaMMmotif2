/*
 * ScoreSeqSet.cpp
 *
 *  Created on: Dec 13, 2016
 *      Author: wanwan
 */

#include "ScoreSeqSet.h"

#include <float.h>		// -FLT_MAX

ScoreSeqSet::ScoreSeqSet( Motif* motif, BackgroundModel* bg, std::vector<Sequence*> seqSet ){

	for( int k = 0; k < Global::Yk; k++ ){
		Y_.push_back( ipow( Alphabet::getSize(), k ) );
	}

	motif_ = motif;
	bg_ = bg;
	seqSet_ = seqSet;

}


ScoreSeqSet::~ScoreSeqSet(){

}

// store the log odds scores at all positions of each sequence
void ScoreSeqSet::score(){

	int K = Global::modelOrder;
	int W = motif_->getW();

	// pre-calculate log odds scores given motif and bg model
	motif_->calculateS( bg_->getV() );
	float** s = motif_->getS();

	mops_scores_.resize( seqSet_.size() );

	for( size_t n = 0; n < seqSet_.size(); n++ ){

		int 	LW1 = seqSet_[n]->getL() - W + 1;
		int* 	kmer = seqSet_[n]->getKmer();
		float 	maxScore = -FLT_MAX;

		for( int i = 0; i < LW1; i++ ){

			float logOdds = 0.0f;

			for( int j = 0; j < W; j++ ){

				int y = ( kmer[i+j] >= 0 ) ? kmer[i+j] % Y_[K+1] : -1;

				logOdds += ( y >= 0 ) ? s[y][j] : 0;

			}

			// take all the log odds scores for MOPS model:
			mops_scores_[n].push_back( logOdds );

			// take the largest log odds score for ZOOPS model:
			maxScore = ( logOdds > maxScore ) ? logOdds : maxScore;
		}
		zoops_scores_.push_back( maxScore );
	}
}

std::vector<std::vector<float>> ScoreSeqSet::getMopsScores(){
	return mops_scores_;
}

std::vector<float> ScoreSeqSet::getZoopsScores(){
	return zoops_scores_;
}

void ScoreSeqSet::write( int N, float cutoff ){

	/**
	 * save log odds scores in one flat file:
	 * posSequenceBasename.logOdds
	 */

	bool 	first_hit = true;
	int 	end; 				// end of motif match

	std::string opath = std::string( Global::outputDirectory ) + '/'
			+ Global::posSequenceBasename +  "_motif_" + std::to_string( N+1 ) + ".logOdds";

	std::ofstream ofile( opath );

	for( size_t n = 0; n < seqSet_.size(); n++ ){
		first_hit = true;
		int seqlen = seqSet_[n]->getL();
		if( !Global::ss ){
			seqlen = seqlen / 2;
		}
		int LW1 = seqSet_[n]->getL() - motif_->getW() + 1;
		for( int i = 0; i < LW1; i++ ){

			if( mops_scores_ [n][i] > cutoff ){
				if( first_hit ){
					// >header:sequence_length
					ofile << '>' << seqSet_[n]->getHeader() <<  ':' << seqlen << std::endl;
					first_hit = false;
				}
				// start:end:score:strand:sequence_matching
				end = i + motif_->getW()-1;

				ofile << i << ':' << end << ':' << std::setprecision( 3 ) << mops_scores_[n][i] << ':' <<
						( ( i < seqlen ) ? '+' : '-' ) << ':' ;
				for( int m = i; m <= end; m++ ){
					ofile << Alphabet::getBase( seqSet_[n]->getSequence()[m] );
				}
				ofile << std::endl;
			}
		}

		//ofile << std::endl;
	}

}
