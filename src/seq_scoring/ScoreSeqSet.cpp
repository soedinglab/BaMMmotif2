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

    size_t seqN = seqSet_.size();
	mops_scores_.resize( seqN );

#pragma omp parallel for
	for( size_t n = 0; n < seqN; n++ ){

		size_t 	LW1 = seqSet_[n]->getL() - W + 1;
		size_t* kmer = seqSet_[n]->getKmer();
        mops_scores_[n].resize( LW1 );    // important for code paralleling
		float 	maxScore = -FLT_MAX;

        size_t z_i = 0;
		for( size_t i = 0; i < LW1; i++ ){
			float logOdds = 0.0f;
			for( size_t j = 0; j < W; j++ ){
				size_t y = kmer[i+j] % Y_[K+1];
				logOdds += s[y][j];
			}
			// take all the log odds scores for MOPS model:
			mops_scores_[n][i] = logOdds;

			// take the largest log odds score for ZOOPS model:
            if( logOdds > maxScore ){
                maxScore = logOdds;
                z_i = i;
            }
        }
		zoops_scores_.push_back( maxScore );
        z_.push_back( z_i );
	}
}

// compute p_values for motif scores based on negative sequence scores
void ScoreSeqSet::calcPvalues( std::vector<std::vector<float>> pos_scores,
                               std::vector<float> neg_all_scores ){

	/**
	 * calculate P-values for motif occurrences
	 */
/*

    size_t posN = seqSet_.size();
    size_t negN = neg_all_scores.size();
    mops_p_values_.resize( posN );
    mops_e_values_.resize( posN );

    float eps = 1.0e-5;

    // sort negative set scores in descending order
    std::sort( neg_all_scores.begin(),
               neg_all_scores.end(),
               std::greater<float>() );

    // get the top n-th score from the negative set
    size_t nTop = std::min( 100, ( int )negN / 10 );
    float S_ntop = neg_all_scores[nTop-1];

    // calculate the rate parameter lambda
    float lambda = 0.f;
	for( size_t n = 0; n < nTop; n++ ){
		lambda += ( neg_all_scores[n] - S_ntop );
	}
	lambda = lambda / ( float )nTop;

#pragma omp parallel for
	for( size_t n = 0; n < seqSet_.size(); n++ ){

		size_t LW1 = seqSet_[n]->getL() - motif_->getW() + 1;

		for( size_t i = 0; i < LW1; i++ ){

            float Sl = pos_scores[n][i];
            // count the accumulated number of scores from the negative set up to rank l
            size_t FPl = std::distance( std::upper_bound( neg_all_scores.begin(),
                                                          neg_all_scores.end(), Sl ),
                                        neg_all_scores.end() );

            float p_value;
			if( FPl == negN ){
                // when Sl is lower than the worst negative score:
				p_value = 1.f;

			} else if( FPl < 10 and fabs( lambda ) > eps ){
			    // when only few or no negatives are higher than S_l:
				p_value = float( nTop ) / ( float )negN * expf( - ( Sl - S_ntop ) / lambda );

                if( p_value >= 1.f ){
                    std::cout << "n=" << n << ",i=" << i
                              << "pval=" << p_value
                              << ", nTop=" << nTop
                              << ", negN=" << negN
                              << ", FPl = " << FPl
                              << ", Sl=" << Sl
                              << ", S_ntop=" << S_ntop
                              << ", lambda=" << lambda
                              << std::endl;
                }

			} else {
				// when Sl_higher and Sl_lower can be defined:
				float SlHigher = neg_all_scores[negN-FPl-1];
				float SlLower = neg_all_scores[negN-FPl];
				p_value = ( ( float )FPl + ( SlHigher - Sl + eps ) / ( SlHigher - SlLower + eps ) )
                          / ( float )negN;

                if( p_value >= 1.f ){
                    std::cout << "n=" << n << ",i=" << i
                              << "pval=" << p_value
                              << ", eps=" << eps
                              << ", negN=" << negN
                              << ", Sl=" << Sl
                              << ", S_higher=" << SlHigher
                              << ", S_lower=" << SlLower
                              << ", FPl=" << FPl
                              << std::endl;
                }
			}

            assert( p_value <= 1.1f );

            mops_p_values_[n].push_back( p_value );
            mops_e_values_[n].push_back( p_value * ( float )posN );
		}
	}
*/

    size_t seqN = seqSet_.size();
    size_t negN = neg_all_scores.size();
    mops_p_values_.resize( seqN );
    mops_e_values_.resize( seqN );
    for( size_t n = 0; n < seqN; n++ ){
        mops_p_values_[n].resize( pos_scores[n].size() );
        mops_e_values_[n].resize( pos_scores[n].size() );
    }

    // create a vector to store the positive sequence lengths
    std::vector<size_t> seql;
    seql.push_back( pos_scores[0].size() );
    for( size_t n = 1; n < seqN; n++ ){
        seql.push_back( seql[n-1] + pos_scores[n].size() );
    }

    // get all scores from positive set
    std::vector<float> pos_all_scores;
    for( size_t n = 0; n < seqN; n++ ){
        for( size_t i = 0; i < pos_scores[n].size(); i++ ){
            pos_all_scores.push_back( pos_scores[n][i] );
        }
    }
    size_t posN = pos_all_scores.size();

    // get the permutation of index after sorting all positive scores in descending order
    std::vector<size_t> pidx_pos = sortIndices( pos_all_scores, true );

    // get the permutation of index after sorting all negative scores in descending order
    std::vector<size_t> pidx_neg = sortIndices( neg_all_scores, true );

    // pre-calculation:
    // get the top n-th score from the negative set
    size_t nTop = std::min( 100, (int)negN / 10 );
    float S_ntop = neg_all_scores[pidx_neg[nTop]];
    // calculate the rate parameter lambda
    float lambda = 0.f;
    for( size_t n = 0; n < nTop; n++ ){
        lambda += ( neg_all_scores[pidx_neg[n]] - S_ntop );
    }
    lambda = lambda / ( float )nTop;
    assert( lambda > 0.0f );

    float  eps = 1.0e-5;   // to avoid 0 in the denominator
    size_t cScore_pos = 0;
    size_t cScore_neg = 0;
    size_t n = 0;
    size_t t = 0;
    for( size_t l = 0; l < posN+negN; l++ ){
        // calculate the accumulated number of scores for both positive and negative score lists
        // Note: for ties (equal scores), the negative scores are always ranked before the positive scores
        if( ( neg_all_scores[pidx_neg[cScore_neg]] >= pos_all_scores[pidx_pos[cScore_pos]]
              and cScore_neg < negN ) or cScore_pos == posN ){
            cScore_neg ++;

        } else {
            // take the accumulated score for negative set as false positive of entry l
            size_t FPl = cScore_neg;

            float Sl = pos_all_scores[pidx_pos[cScore_pos]];
            float Sl_higher = neg_all_scores[pidx_neg[cScore_neg]];
            float Sl_lower  = neg_all_scores[pidx_neg[cScore_neg]];

            size_t i = 1;
            while( Sl_higher < Sl ){
                if( cScore_neg <= i ){ break; }
                Sl_higher = neg_all_scores[pidx_neg[cScore_neg-i]];
                i++;
            }

            size_t j = 1;
            while( Sl_lower > Sl ){
                if( cScore_neg+j > negN ){ break; }
                Sl_lower = neg_all_scores[pidx_neg[cScore_neg+j]];
                j++;
            }

            // calculate p-values for entry l
            float pVal;
            if( FPl > 10 or fabs( lambda ) < eps ) {
                pVal = ((float) FPl + (Sl_higher - Sl) / (Sl_higher - Sl_lower + eps)) / (float) negN;
            } else {
                pVal = nTop * expf( ( S_ntop - Sl ) / lambda ) / (float)negN;
            }

            //assert( pVal <= 1.00001f and pVal >= 0.f );
            n = std::distance( seql.begin(),
                               std::lower_bound(seql.begin(), seql.end(), pidx_pos[cScore_pos]) );

            if( pidx_pos[cScore_pos] == seql[n] and n < seql[seqN-1] ){
                n++;
            }
            t = ( n == 0) ? pidx_pos[cScore_pos] : pidx_pos[cScore_pos] - seql[n-1];

            mops_p_values_[n][t] = pVal;
            mops_e_values_[n][t] = pVal * posN;

            // add one to the accumulated score of positive set
            cScore_pos ++;
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

void ScoreSeqSet::printLogOdds(){
    for( size_t n = 0; n < seqSet_.size(); n++ ){
        std::cout << "seq " << n << ":" << std::endl;
        std::cout << zoops_scores_[n] << '\t';
        for( size_t i = 0; i < mops_scores_[n].size(); i++ ){
            std::cout << mops_scores_[n][i] << '\t';
        }
        std::cout << std::endl;
    }
}

void ScoreSeqSet::write( char* odir, std::string basename, float pvalCutoff, bool ss ){

    /**
	 * save log odds scores in one flat file:
	 * basename.occurrence
	 */

    assert( pval_is_calulated_ );

    // fix eValCutoff for small sequence sets
//    float eValCutoff = 0.1f;

	size_t 	end; 				// end of motif match

	std::string opath = std::string( odir )  + '/' + basename + ".occurrence";

	std::ofstream ofile( opath );

    // add a header to the results
    ofile << "seq\tlength\tstrand\tstart..end\tpattern\tp-value\te-value\tlogOddsScore" << std::endl;

	for( size_t n = 0; n < seqSet_.size(); n++ ){
		size_t seqlen = seqSet_[n]->getL();
		if( !ss ){
			seqlen = ( seqlen - 1 ) / 2;
		}
		size_t LW1 = seqSet_[n]->getL() - motif_->getW() + 1;

        for( size_t i = 0; i < LW1; i++ ){

//            if( mops_p_values_[n][i] < pvalCutoff and mops_scores_[n][i] >= 13.f ){
            if( mops_p_values_[n][i] < pvalCutoff ){
//			  if( mops_e_values_[n][i] < eValCutoff ){
                // >header:sequence_length
                ofile << seqSet_[n]->getHeader() << '\t' << seqlen << '\t';

				// start:end:score:strand:sequence_matching
				end = i + motif_->getW();

				ofile << ( ( i < seqlen ) ? '+' : '-' ) << '\t'
                      << i+1 << ".." << end << '\t';
				for( size_t m = i; m < end; m++ ){
					ofile << Alphabet::getBase( seqSet_[n]->getSequence()[m] );
				}
				ofile << '\t' << std::setprecision( 3 )
                      << mops_p_values_[n][i] << '\t'
                      << mops_e_values_[n][i] << '\t'
                      << mops_scores_[n][i] << std::endl;
			}
		}
	}

}

void ScoreSeqSet::writeLogOdds( char* odir, std::string basename, bool ss ){

    /**
	 * save log odds scores in one flat file:
	 * basename.logOddsZoops
	 */

    size_t 	end; 				// end of motif match

    std::string opath = std::string( odir )  + '/' + basename + ".logOddsZoops";

    std::ofstream ofile( opath );

    // add a header to the results
    ofile << "seq\tlength\tstrand\tstart..end\tpattern\tzoops_score" << std::endl;

    for( size_t n = 0; n < seqSet_.size(); n++ ){
        size_t seqlen = seqSet_[n]->getL();
        if( !ss ){
            seqlen = ( seqlen - 1 ) / 2;
        }

        // >header:sequence_length
        ofile << seqSet_[n]->getHeader() << '\t' << seqlen << '\t';

        // start:end:score:strand:sequence_matching
        end = z_[n] + motif_->getW();

        ofile << ( ( z_[n] < seqlen ) ? '+' : '-' ) << '\t'
              << z_[n]+1 << ".." << end << '\t';
        for( size_t m = z_[n]; m < end; m++ ){
            ofile << Alphabet::getBase( seqSet_[n]->getSequence()[m] );
        }
        ofile << '\t' << std::setprecision( 3 ) << zoops_scores_[n] << std::endl;


    }

}
