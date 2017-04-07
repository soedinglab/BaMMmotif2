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

	s_ = ( float** )calloc( Y_[Global::modelOrder+1], sizeof( float* ) );
	for( int y = 0; y < Y_[Global::modelOrder+1]; y++ ){
		s_[y] = ( float* )calloc( motif_->getW(), sizeof( float ) );
	}
	seqSet_ = seqSet;

	// allocate memory for scores_[n][i]
	scores_.resize( seqSet_.size() );

}

ScoreSeqSet::~ScoreSeqSet(){
	for( int y = 0; y < Y_[Global::modelOrder+1]; y++ ){
		free( s_[y]);
	}
	free( s_ );
}

// store the log odds scores at all positions of each sequence
void ScoreSeqSet::score(){

	int 	K = Global::modelOrder;
	int 	K_bg = ( Global::bgModelOrder < K ) ? Global::bgModelOrder : K;
	int 	W = motif_->getW();

	//	float 	maxScore;								// maximal logOddsScore over all positions for each sequence

	for( int y = 0; y < Y_[Global::modelOrder+1]; y++ ){
	    int y_bg = y % Y_[K_bg+1];
		for( int j = 0; j < W; j++ ){
			s_[y][j] = logf( motif_->getV()[K][y][j] ) - logf( bg_->getV()[K_bg][y_bg] );
		}
	}

	for( size_t n = 0; n < seqSet_.size(); n++ ){
		int LW1 = seqSet_[n]->getL() - W + 1;
		int* kmer = seqSet_[n]->getKmer();

		float logOdds;

		for( int i = 0; i < LW1; i++ ){
			logOdds = 0.0f;

			for( int j = 0; j < W; j++ ){
				int y = ( kmer[i+j] >= 0 ) ? kmer[i+j] % Y_[K+1] : -1;
				logOdds += ( y >= 0 ) ? s_[y][j] : 0;
			}

			// take all the log odds scores for MOPS model:
			scores_[n].push_back( logOdds );

			// take the largest log odds score for ZOOPS model:
			//if( logOdds[i] > maxScore ){
			//	maxScore = logOdds[i];
			//}
		}
		//scores_[1].push_back( maxScore );
	}
}

std::vector<std::vector<float>> ScoreSeqSet::getScores(){
	return scores_;
}

void ScoreSeqSet::write( int N, float cutoff ){

	/**
	 * save log odds scores in one flat file:
	 * posSequenceBasename.logOdds
	 */

	bool first_hit = true;
	int end; 				// end of motif match

	std::string opath = std::string( Global::outputDirectory ) + '/'
			+ Global::posSequenceBasename +  "_motif_" + std::to_string( N+1 ) + ".logOdds";

	std::ofstream ofile( opath );

	for( size_t n = 0; n < seqSet_.size(); n++ ){
		first_hit = true;
		int seqlen = seqSet_[n]->getL();
		if( !Global::ss ){
			seqlen = seqlen/2;
		}
		int LW1 = seqSet_[n]->getL() - motif_->getW() + 1;
		for( int i = 0; i < LW1; i++ ){

			if( scores_ [n][i] > cutoff ){
				if ( first_hit ){
					// >header:sequence_length
					ofile << '>' << seqSet_[n]->getHeader() <<  ':' << seqlen << std::endl;
					first_hit = false;
				}
				// start:end:score:strand:sequence_matching
				end = i + motif_->getW()-1;

				ofile << i << ':' << end << ':' << std::setprecision( 3 ) << scores_[n][i] << ':' <<
						((i < seqlen) ? '+' : '-') << ':' ;
				for( int m = i; m <= end; m++ ){
					ofile << Alphabet::getBase( seqSet_[n]->getSequence()[m] );
				}
				ofile << std::endl;
			}
		}

		//ofile << std::endl;
	}

}
