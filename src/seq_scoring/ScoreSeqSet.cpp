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

    // sort negative set scores in ascending order
    std::sort( neg_all_scores.begin(), neg_all_scores.end(), std::less<float>() );

    // get the top n-th score from the negative set
    size_t nTop = std::min( 100, ( int )negN / 10 );
    float S_ntop = neg_all_scores[nTop];

    // calculate the rate parameter lambda
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
			} else if( FPl < 10 and fabs( lambda ) > eps ){
			    // when only few or no negatives are higher than S_l:
				p_value = float( nTop ) / ( float )negN * expf( - ( Sl - S_ntop ) / lambda );

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

    /*
    size_t posN = seqSet_.size();
    size_t negN = neg_all_scores.size();
    mops_p_values_.resize( posN );
    mops_e_values_.resize( posN );

    float eps = 1.0e-5;

    std::vector<float> pos_all_scores;
    for( size_t n = 0; n < pos_scores.size(); n++ ){
        for( size_t i = 0; i < pos_scores[n].size(); i++ ){
            pos_all_scores.push_back( pos_scores[n][i] );
        }
    }

    std::vector<float> all_scores;
    all_scores = pos_all_scores;
    for( size_t m = 0; m < neg_all_scores.size(); m++ ){
        all_scores.push_back( neg_all_scores[m] );
    }

    // get the permutation of index after sorting all positive and negative scores
    // jointly in descending order
    std::vector<size_t> pidx_all = sortIndices( all_scores );

    // get the permutation of index after sorting all positive scores in descending order
    std::vector<size_t> pidx_pos = sortIndices( pos_all_scores );

    // get the permutation of index after sorting all negative scores in descending order
    std::vector<size_t> pidx_neg = sortIndices( neg_all_scores );

    size_t cScore_pos = 0;
    size_t cScore_neg = 0;

    std::vector<size_t> FP;
    std::vector<float> pValues;
    std::vector<float> eValues;
    // pre-calculation:
    // get the top n-th score from the negative set
    size_t nTop = std::min( 100, (int)negN / 10 );
    float S_ntop = neg_all_scores[pidx_neg[nTop-1]];
    // calculate the rate parameter lambda
    float lambda = 0.f;
    for( size_t n = 0; n < nTop; n++ ){
        lambda += ( neg_all_scores[pidx_neg[n]] - S_ntop );
    }
    lambda = lambda / ( float )nTop;

    std::cout << S_ntop << '\t' << lambda << std::endl;
    std::cout << pos_all_scores.size() << '\t' << neg_all_scores.size() << std::endl;

    for( size_t l = 0; l < all_scores.size(); l++ ){
        // calculate the accumulated number of scores for both positive and negative score lists
        // Note: for ties (equal scores), the negative scores are always ranked before the positive scores
        if( ( neg_all_scores[pidx_neg[cScore_neg]] >= pos_all_scores[pidx_pos[cScore_pos]] and cScore_neg < negN )
            or cScore_pos == pos_all_scores.size() ){
            cScore_neg ++;
        } else {
            cScore_pos ++;
        }

        // take the accumulated score for negative set as false positive of entry l
        FP.push_back( cScore_neg );

        float Sl = all_scores[pidx_all[l]];
        float Sl_higher = neg_all_scores[pidx_neg[cScore_neg] - 1];
        float Sl_lower = neg_all_scores[pidx_neg[cScore_neg]];

        // calculate p-values for entry l
        float pVal;
        if( FP[l] > 10 or fabs( lambda ) < eps ) {
            pVal = ((float) FP[l] + (Sl_higher - Sl) / (Sl_higher - Sl_lower + eps)) / (float) negN;
        } else {
            pVal = nTop * expf( ( S_ntop - Sl ) / lambda ) / (float)negN;
        }
        pValues.push_back(pVal);

        // if(l < 100) std::cout << pVal << std::endl;
        //std::cout << l << '\t';
    }

    // assign p-value to each position on the sequence
    size_t site = 0;

    for( size_t n = 0; n < posN; n++ ){

		size_t LW1 = seqSet_[n]->getL() - motif_->getW() + 1;

		for( size_t i = 0; i < LW1; i++ ){
            mops_p_values_[n].push_back( pValues[pidx_all[site]] );
            mops_e_values_[n].push_back( pValues[pidx_all[site]] * posN );
            site ++;
		}
	}
*/
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
	 * basename.occurrence
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
