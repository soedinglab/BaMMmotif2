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
    pval_is_calulated_ = false;
}

ScoreSeqSet::~ScoreSeqSet(){

}


void ScoreSeqSet::calcLogOdds(){

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

// compute p_values for motif scores based on negative sequence scores
void ScoreSeqSet::calcPvalues( std::vector<std::vector<float>> pos_scores, std::vector<float> neg_all_scores ){

	/**
	 * calculate P-values for motif occurrences
	 */
    size_t posN = seqSet_.size();
    size_t negN = neg_all_scores.size();
    mops_p_values_.resize( posN );
    mops_e_values_.resize( posN );

    float eps = 1.0e-5;
    size_t nTop = std::min( 100, ( int )negN / 10 );

    // sort negative set scores in ascending order
    std::sort( neg_all_scores.begin(), neg_all_scores.end(), std::less<float>() );

    float S_ntop = neg_all_scores[nTop];
    float lambda = 0.f;
	for( size_t n = 0; n < nTop; n++ ){
		lambda += ( neg_all_scores[n] - S_ntop );
	}
	lambda = lambda / ( float )nTop;

	for( size_t n = 0; n < seqSet_.size(); n++ ){

		size_t LW1 = seqSet_[n]->getL() - motif_->getW() + 1;

		for( size_t i = 0; i < LW1; i++ ){

            float Sl = pos_scores[n][i];
            // count the accumulated number of scores from the negative set up to rank l
            size_t FPl = std::distance( std::upper_bound( neg_all_scores.begin(), neg_all_scores.end(), Sl ),
                                        neg_all_scores.end() );

            float p_value;
			if( FPl == negN ){
                // when Sl is lower than the worst negative score:
				p_value = 1.f;
			} else if( FPl < 10 ){
                // when only few or no negative scores are higher than Sl:
				p_value = float( nTop ) / ( float )negN * exp( - ( Sl - S_ntop ) / lambda );
			} else {
				// when Sl_higher and Sl_lower can be defined:
				float SlHigher = neg_all_scores[negN-FPl-1];
				float SlLower = neg_all_scores[negN-FPl];
				p_value = ( ( float )FPl + ( SlHigher - Sl + eps ) / ( SlHigher - SlLower + eps ) ) / ( float )negN;
			}
            mops_p_values_[n].push_back( p_value );
            mops_e_values_[n].push_back( p_value * ( float )posN );
		}
	}

    pval_is_calulated_ = true;
}

std::vector<std::vector<float>> ScoreSeqSet::getMopsScores(){
	return mops_scores_;
}

std::vector<float> ScoreSeqSet::getZoopsScores(){
	return zoops_scores_;
}

void ScoreSeqSet::write( char* odir, std::string basename, float pvalCutoff, bool ss ){
	/**
	 * save log odds scores in one flat file:
	 * posSequenceBasename_motif_N.occurrence
	 */

    assert( pval_is_calulated_ );

	size_t 	end; 				// end of motif match

	std::string opath = std::string( odir )  + '/' + basename + ".occurrence";

	std::ofstream ofile( opath );

	for( size_t n = 0; n < seqSet_.size(); n++ ){
		size_t seqlen = seqSet_[n]->getL();
		if( !ss ){
			seqlen = ( seqlen - 1 ) / 2;
		}
		size_t LW1 = seqSet_[n]->getL() - motif_->getW() + 1;

        for( size_t i = 0; i < LW1; i++ ){

			if( mops_p_values_[n][i] < pvalCutoff ){
                // >header:sequence_length
                ofile << '>' << seqSet_[n]->getHeader() << '\t' << seqlen << '\t';

				// start:end:score:strand:sequence_matching
				end = i + motif_->getW()-1;

				ofile << ( ( i < seqlen ) ? '+' : '-' ) << '\t'
                      << i << ".." << end << '\t';
				for( size_t m = i; m <= end; m++ ){
					ofile << Alphabet::getBase( seqSet_[n]->getSequence()[m] );
				}
				ofile << '\t' << std::setprecision( 3 )
                      << mops_p_values_[n][i] << '\t'
                      << mops_e_values_[n][i] << std::endl;
			}
		}
	}

}
