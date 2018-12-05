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

    size_t seqN = seqSet_.size();
    mops_p_values_.resize( seqN );
    mops_e_values_.resize( seqN );
    mops_scores_.resize( seqN );
    seql_.resize( seqN );
    for( size_t n = 0; n < seqN; n++ ){
        size_t 	LW1 = seqSet_[n]->getL() - motif->getW() + 1;
        mops_scores_[n].resize( LW1 );    // important for code paralleling
        mops_p_values_[n].resize( LW1 );
        mops_e_values_[n].resize( LW1 );
        if(n == 0 ){
            seql_[n] = LW1;
        } else {
            seql_[n] = seql_[n-1] + LW1;
        }
    }

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

    // todo: the parallelization does not work properly given PWM and applied EM
//#pragma omp parallel for
	for( size_t n = 0; n < seqN; n++ ){

		size_t 	LW1 = seqSet_[n]->getL() - W + 1;
		size_t* kmer = seqSet_[n]->getKmer();

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

    size_t seqN = seqSet_.size();

    // copy all scores from positive set into one vector
    std::vector<float> pos_all_scores;
    for (size_t i = 0; i < seqN; i++) {
        pos_all_scores.insert(pos_all_scores.end(),
                              pos_scores[i].begin(),
                              pos_scores[i].end());
    }

    size_t posN = pos_all_scores.size();
    size_t negN = neg_all_scores.size();

    std::cout << "There are " << posN << " positions in the given set." << std::endl;

    float eps = 1.0e-5f; // to avoid 0 in the denominator

    bool run_slow = true;

    if( !run_slow ) {

        // sort negative set scores in descending order
        std::sort(neg_all_scores.begin(),
                  neg_all_scores.end(),
                  std::greater<float>());

        // get the top n-th score from the negative set
        size_t nTop = static_cast<size_t>( std::min(100, (int) negN / 10) );
        float S_ntop = neg_all_scores[nTop-1];

        // calculate the rate parameter lambda
        float lambda = 1.e-5f;  // to avoid lambda from being 0
        for (size_t n = 0; n < nTop; n++) {
            lambda += (neg_all_scores[n] - S_ntop);
        }
        lambda = lambda / (float) nTop;
        assert(lambda > 0.0f);

#pragma omp parallel for
        for (size_t n = 0; n < seqN; n++) {

            size_t LW1 = pos_scores[n].size();

            for (size_t i = 0; i < LW1; i++) {

                float Sl = pos_scores[n][i];
                // count the accumulated number of scores from the negative set up to rank l
                size_t FPl = static_cast<size_t>( std::distance(
                        std::upper_bound( neg_all_scores.begin(),
                                          neg_all_scores.end(), Sl),
                        neg_all_scores.end()));

                float pVal;

                if (FPl == negN) {
                    // when Sl is lower than the worst negative score:
                    pVal = 1.f;

                } else if (FPl < 10 and fabs(lambda) > eps) {
                    // when only few or no negatives are higher than S_l:
                    pVal = float(nTop) * expf((S_ntop - Sl) / lambda)
                           / (float) negN;

                    if (pVal >= 1.f) {
                        std::cout << "n=" << n << ",i=" << i
                                  << ", pval=" << pVal
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
                    float Sl_higher = neg_all_scores[negN-FPl];
                    float Sl_lower = neg_all_scores[negN-FPl];

                    size_t ii = 1;
                    while (Sl_higher < Sl) {
                        if (Sl_higher <= ii) { break; }
                        Sl_higher = neg_all_scores[negN-FPl-ii];
                        ii++;
                    }

                    size_t j = 1;
                    while (Sl_lower > Sl) {
                        if (Sl_lower + j > negN) { break; }
                        Sl_lower = neg_all_scores[negN-FPl+j];
                        j++;
                    }
                    pVal = ((float) FPl + (Sl_higher - Sl + eps) / (Sl_higher - Sl_lower + eps))
                              / (float) negN;

                    if (pVal >= 1.f) {
                        std::cout << "n=" << n << ",i=" << i
                                  << ", pval=" << pVal
                                  << ", eps=" << eps
                                  << ", negN=" << negN
                                  << ", Sl=" << Sl
                                  << ", S_higher=" << Sl_higher
                                  << ", S_lower=" << Sl_lower
                                  << ", FPl=" << FPl
                                  << std::endl;
                    }
                }

                //assert(pVal <= 1.1f);

                mops_p_values_[n][i] = pVal;
                mops_e_values_[n][i] = pVal * posN;
            }
        }
    } else {

        // get the permutation of index after sorting all positive scores
        // in descending order
        std::vector<size_t> pidx_pos = sortIndices(pos_all_scores, true);

        // get the permutation of index after sorting all negative scores
        // in descending order
        std::vector<size_t> pidx_neg = sortIndices(neg_all_scores, true);

        // pre-calculation:
        // get the top n-th score from the negative set
        size_t nTop = static_cast<size_t>(std::min(100, (int) negN / 10));
        float S_ntop = neg_all_scores[pidx_neg[nTop]];
        // calculate the rate parameter lambda
        float lambda = 1.e-5f;  // to avoid lambda from being 0
        for (size_t c = 0; c < nTop; c++) {
            lambda += (neg_all_scores[pidx_neg[c]] - S_ntop);
        }
        lambda = lambda / (float) nTop;
        assert(lambda > 0.0f);

        size_t cScore_pos = 0;
        size_t cScore_neg = 0;
        size_t n = 0;
        size_t t = 0;

        for (size_t l = 0; l < posN + negN; l++) {
            // calculate the accumulated number of scores for both
            // positive and negative score lists
            // Note: for ties (equal scores), the negative scores are
            // always ranked before the positive scores
            if( ( neg_all_scores[pidx_neg[cScore_neg]] >=
                          pos_all_scores[pidx_pos[cScore_pos]]
                  or cScore_pos == posN-1 ) and cScore_neg < negN ) {

                cScore_neg++;

//                std::cout << l << '\t'
//                          << pos_all_scores[pidx_pos[cScore_pos]] << '\t'
//                          << neg_all_scores[pidx_neg[cScore_neg]] << '\t'
//                          << cScore_pos << '\t'
//                          << cScore_neg << '\n';
            } else {
                // take the accumulated score for negative set as
                // false positive of entry l
                size_t FPl = cScore_neg;

//                std::cout << l << '\t'
//                          << pos_all_scores[pidx_pos[cScore_pos]] << '\t'
//                          << neg_all_scores[pidx_neg[cScore_neg]] << '\t'
//                          << cScore_pos << '\t'
//                          << cScore_neg << " (+)" <<'\n';

                float Sl = pos_all_scores[pidx_pos[cScore_pos]];

                // calculate p-values for entry l
                float pVal;

                if ( FPl == negN ){
                    pVal = 1.f;
                } else if(FPl > 10 or fabs(lambda) < eps) {

                    float Sl_higher = neg_all_scores[pidx_neg[cScore_neg]];
                    float Sl_lower = neg_all_scores[pidx_neg[cScore_neg]];

                    size_t i = 1;
                    while (Sl_higher < Sl) {
                        if (cScore_neg <= i) { break; }
                        Sl_higher = neg_all_scores[pidx_neg[cScore_neg-i]];
                        i++;
                    }

                    size_t j = 1;
                    while (Sl_lower > Sl) {
                        if (cScore_neg + j > negN) { break; }
                        Sl_lower = neg_all_scores[pidx_neg[cScore_neg+j]];
                        j++;
                    }

                    pVal = ((float) FPl + (Sl_higher - Sl) / (Sl_higher - Sl_lower + eps))
                           / (float) negN;
                } else {
                    pVal = nTop * expf((S_ntop - Sl) / lambda) / (float) negN;
                }

                assert( pVal <= 1.00001f and pVal >= 0.f );

                n = static_cast<size_t>( std::distance(seql_.begin(),
                                                       std::lower_bound(seql_.begin(),
                                                                        seql_.end(),
                                                                        pidx_pos[cScore_pos])) );

                if (pidx_pos[cScore_pos] == seql_[n] and n < seql_[seqN-1]) {
                    n++;
                }

                if(n == 0){
                    t = pidx_pos[cScore_pos];
                } else {
                    t = pidx_pos[cScore_pos] - seql_[n-1];
                }

                mops_p_values_[n][t] = pVal;
                mops_e_values_[n][t] = pVal * posN;

                // add one to the accumulated score of positive set
                cScore_pos++;
            }
        }

    }

    pval_is_calulated_ = true;

//    // clear the content
//    pidx_pos.clear();
//    pidx_neg.clear();
//    pos_all_scores.clear();

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

//            if( mops_scores_[n][i] >= 0.f ){
//            if( mops_scores_[n][i] >= 10.f ){
            if( mops_p_values_[n][i] < pvalCutoff ){
//			  if( mops_e_values_[n][i] < 0.1f ){
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
