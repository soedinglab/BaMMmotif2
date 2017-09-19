/*
 * ScoreSeqSet.cpp
 *
 *  Created on: Dec 13, 2016
 *      Author: wanwan
 */

#include "ScoreSeqSet.h"

#include <float.h>		// -FLT_MAX

ScoreSeqSet::ScoreSeqSet( Motif* motif, BackgroundModel* bg, std::vector<Sequence*> seqSet ){

	motif_	= motif;
	bg_ 	= bg;
	seqSet_	= seqSet;
	Y_ 		= motif->getY();

}

ScoreSeqSet::~ScoreSeqSet(){

}


void ScoreSeqSet::score(){

	/**
	 * store the log odds scores at all positions of each sequence
	 */

	size_t K = motif_->getK();
	size_t W = motif_->getW();
	size_t K_bg = ( bg_->getOrder() < K ) ? bg_->getOrder() : K;
	// pre-calculate log odds scores given motif and bg model
	motif_->calculateLogS( bg_->getV(), K_bg );
	float** s = motif_->getS();

	mops_scores_.resize( seqSet_.size() );

	for( size_t n = 0; n < seqSet_.size(); n++ ){

		size_t 	LW1 = seqSet_[n]->getL() - W + 1;
		size_t* kmer = seqSet_[n]->getKmer();
		float 	maxScore = -FLT_MAX;

		for( size_t i = 0; i < LW1; i++ ){

			float logOdds = 0.0f;

			for( size_t j = 0; j < W; j++ ){

				size_t y = kmer[i+j] % Y_[K+1];

				logOdds += s[y][j];

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

void ScoreSeqSet::write( char* odir, std::string basename, float cutoff, bool ss ){
	/**
	 * save log odds scores in one flat file:
	 * posSequenceBasename_motif_N.occurrence
	 */

	bool 	first_hit = true;
	size_t 	end; 				// end of motif match

	std::string opath = std::string( odir )  + '/' + basename + ".occurrence";

	std::ofstream ofile( opath );

	for( size_t n = 0; n < seqSet_.size(); n++ ){
		first_hit = true;
		size_t seqlen = seqSet_[n]->getL();
		if( !ss ){
			seqlen = ( seqlen - 1 ) / 2;
		}
		size_t LW1 = seqSet_[n]->getL() - motif_->getW() + 1;
		for( size_t i = 0; i < LW1; i++ ){

			if( mops_scores_ [n][i] > cutoff ){
				if( first_hit ){
					// >header:sequence_length
					ofile << '>' << seqSet_[n]->getHeader() << ':' << seqlen << std::endl;
					first_hit = false;
				}
				// start:end:score:strand:sequence_matching
				end = i + motif_->getW()-1;

				ofile << i << ".." << end << ':' << std::setprecision( 3 )
						<< mops_scores_[n][i] << ':' <<
						( ( i < seqlen ) ? '+' : '-' ) << ':' ;
				for( size_t m = i; m <= end; m++ ){
					ofile << Alphabet::getBase( seqSet_[n]->getSequence()[m] );
				}
				ofile << std::endl;
			}
		}
	}

}
