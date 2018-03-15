#include "FDR.h"

FDR::FDR( std::vector<Sequence*> posSeqs, std::vector<Sequence*> negSeqs, float q,
          Motif* motif, BackgroundModel* bgModel, size_t cvFold,
          bool mops, bool zoops, bool savePRs,
          bool savePvalues, bool saveLogOdds){

	posSeqs_	= posSeqs;
	negSeqs_	= negSeqs;
	q_ 			= q;
	motif_ 		= motif;
    bgModel_    = bgModel;
	cvFold_		= cvFold;
    mops_       = mops;
    zoops_      = zoops;
    savePRs_    = savePRs;
    savePvalues_= savePvalues;
    saveLogOdds_= saveLogOdds;
	occ_frac_	= 0.0f;
	occ_mult_	= 0.0f;

}

FDR::~FDR(){

}

void FDR::evaluateMotif( bool EMoptimize, bool CGSoptimize, bool optimizeQ, bool advanceEM, float frac ){

	std::vector<std::vector<float>> mops_scores;
	std::vector<float> 				zoops_scores;
    float updatedQ = q_;  // obtain the updated q

    /**
	 * Cross validation
	 */
	for( size_t fold = 0; fold < cvFold_; fold++ ){

		// deep copy the initial motif
		Motif* motif = new Motif( *motif_ );

		/**
		 * Draw sequences for each training, test and negative sets
		 */
		std::vector<Sequence*> testSet;
		std::vector<Sequence*> trainSet;
        std::vector<Sequence*> negSet;
		for( size_t n = 0; n <= posSeqs_.size()-cvFold_; n+=cvFold_ ){
			for( size_t f = 0; f < cvFold_; f++ ){
				if( f != fold ){
					trainSet.push_back( posSeqs_[n+f] );
				} else {
					testSet.push_back( posSeqs_[n+f] );
				}
			}
		}
		for( size_t n = 0; n <= negSeqs_.size()-cvFold_; n+=cvFold_ ){
            negSet.push_back( negSeqs_[n] );
		}

		/**
		 * Training
		 */
		// learn motif from each training set
		if( EMoptimize ){
			EM model( motif, bgModel_, trainSet, q_, optimizeQ, false, frac );
            if( advanceEM ){
                model.advance();
            } else {
                model.optimize();
            }
            updatedQ = model.getQ();
		} else if ( CGSoptimize ){
			GibbsSampling model( motif, bgModel_, trainSet, q_, optimizeQ );
			model.optimize();
            updatedQ = model.getQ();
		}

		/**
		 * Testing
		 */
		// score positive test sequences with (learned) motif
		ScoreSeqSet score_testset( motif, bgModel_, testSet );
        score_testset.calcLogOdds();

		if( mops_ ){
			mops_scores = score_testset.getMopsScores();
			for( size_t n = 0; n < testSet.size(); n++ ){
				posScoreAll_.insert( std::end( posScoreAll_ ),
                                     std::begin( mops_scores[n] ),
                                     std::end( mops_scores[n] ) );
			}
		}

		if( zoops_ ){
			zoops_scores = score_testset.getZoopsScores();
			posScoreMax_.insert( std::end( posScoreMax_ ),
                                 std::begin( zoops_scores ),
                                 std::end( zoops_scores ) );
		}

		// score negative sequence set
		ScoreSeqSet score_negset( motif, bgModel_, negSet );
        score_negset.calcLogOdds();
		if( mops_ ){
			mops_scores = score_negset.getMopsScores();
			for( size_t n = 0; n < negSet.size(); n++ ){
				negScoreAll_.insert( std::end( negScoreAll_ ),
                                     std::begin( mops_scores[n] ),
                                     std::end( mops_scores[n] ) );
			}
		}
		if( zoops_ ){
			zoops_scores = score_negset.getZoopsScores();
			negScoreMax_.insert( std::end( negScoreMax_ ),
                                 std::begin( zoops_scores ),
                                 std::end( zoops_scores ) );
		}

		if( motif ) 				delete motif;

	}

    // update Q
    q_ = updatedQ;

    // calculate precision and recall
    calculatePR();

    if( savePvalues_ ){

		fprintf( stderr, " ______________________\n"
						"|                      |\n"
						"|  calculate P-values  |\n"
						"|______________________|\n\n" );

		calculatePvalues();
	}

}

void FDR::calculatePR(){

	size_t posN = posSeqs_.size();
	size_t negN = negSeqs_.size();
	float mFold = ( float )negN / ( float )posN;

	// for MOPS model:
	if( mops_ ){
		// Sort log odds scores in descending order
		std::sort( posScoreAll_.begin(), posScoreAll_.end(), std::greater<float>() );
		std::sort( negScoreAll_.begin(), negScoreAll_.end(), std::greater<float>() );

		// Rank and score these log odds score values
		size_t idx_posAll = 0;
		size_t idx_negAll = 0;		// index for arrays storing log odds scores
		float E_TP_MOPS = 0.0f;		// maximal distance between TFP and
									// FP = expectation value of TP
		size_t idx_max = posN+negN;	// index when the distance between TFP and
									// FP reaches maximum; otherwise, set the
									// the number as initial cutoff

		size_t len_all = posScoreAll_.size() + negScoreAll_.size();

		for( size_t i = 0; i < len_all; i++ ){
			if( posScoreAll_[idx_posAll] > negScoreAll_[idx_negAll] || idx_negAll == len_all ){
				idx_posAll++;
			} else {
				idx_negAll++;
			}

			MOPS_TP_.push_back( ( float )idx_posAll - ( float )idx_negAll / mFold );
			MOPS_FP_.push_back( ( float )idx_negAll / mFold );

			if( E_TP_MOPS == MOPS_TP_[i] ){
				idx_max = i;
//				break;				// stop when the distance between TFP and
									// FP reaches maximum
			}
			if( E_TP_MOPS < MOPS_TP_[i] ){
				E_TP_MOPS = MOPS_TP_[i];
			}
		}

		for( size_t i = 0; i < idx_max; i++ ){
			MOPS_FDR_.push_back( MOPS_FP_[i] / ( MOPS_TP_[i] + MOPS_FP_[i] ) );
			MOPS_Rec_.push_back( MOPS_TP_[i] / E_TP_MOPS );
		}

		// the number of motif occurrences per sequence
		occ_mult_ = E_TP_MOPS / ( float )posN;
	}

	// for ZOOPS model:
	if( zoops_ ){
		PN_Pvalue_.clear();

		// Sort log odds scores in descending order
		std::sort( posScoreMax_.begin(), posScoreMax_.end(), std::greater<float>() );
		std::sort( negScoreMax_.begin(), negScoreMax_.end(), std::greater<float>() );

		// Rank and score these log odds score values
		size_t idx_posMax = 0;
		size_t idx_negMax = 0;
		size_t min_idx_pos = 0;
		size_t posN_est = static_cast<size_t>( q_ * ( float )posN );

        // set limit for using the exponential extrapolation for p-value calculation
        size_t n_top = 100 < (negN / 10) ? 100 : (negN / 10);

        float lambda = 0.f;
        for( size_t l = 0; l < n_top; l++ ){
            lambda += negScoreMax_[l] - negScoreMax_[n_top];
        }
        lambda /= n_top;

        float Sl = 0.f;

		for( size_t i = 0; i < posN + negN; i++ ){

            if( idx_posMax <= posN-cvFold_ and idx_negMax <= negN-cvFold_ ){
                if( posScoreMax_[idx_posMax] > negScoreMax_[idx_negMax] or idx_negMax == posN+negN-1){
                    Sl = posScoreMax_[idx_posMax];
                    idx_posMax++;
                } else {
                    Sl = negScoreMax_[idx_negMax];
                    idx_negMax++;
                }
            }

/*
            // stops when TP = FP
            if( ( float )idx_posMax <= ( float )idx_negMax / mFold and i > posN ){
                break;
            }
*/
            float TP = ( float )idx_posMax;
            float FP = ( float )idx_negMax / mFold;
			ZOOPS_TP_.push_back( TP );
			ZOOPS_FP_.push_back( FP );

            float p_value;

            // calculate p-values in two different ways:
            if( Sl < negScoreMax_[n_top] ){
                p_value = ( ( float )idx_negMax + 0.5f ) / ( ( float )negN + 1.0f );

            } else {
                // p-value is calculated by relying on a parametric fit of the exponentially cumulative distribution
                p_value = n_top * expf( ( negScoreMax_[n_top] - Sl ) / lambda ) / negN;
            }

			PN_Pvalue_.push_back( p_value );

			// take the faction of q sequences as real positives
			if( idx_posMax == posN_est ){
				min_idx_pos = i;
			}

			ZOOPS_FDR_.push_back( FP / ( TP + FP ) );
			ZOOPS_Rec_.push_back( TP / ( float )posN );

		}

        // the fraction of motif occurrence
		occ_frac_ = 1.0f - ZOOPS_FP_[min_idx_pos] / ( float )posN;
	}
}

void FDR::calculatePvalues(){

	/**
	 *  calculate p-values for the log odds scores from positive sequence set
	 * 	later used for ranking models using fdrtool R package
	 */

    // for MOPS model:
    if( mops_ ){
        // Sort log odds scores in ascending order
        std::sort( negScoreAll_.begin(), negScoreAll_.end(), std::less<float>() );
        std::sort( posScoreAll_.begin(), posScoreAll_.end(), std::less<float>() );
        for( size_t i = 0; i < posScoreAll_.size(); i++ ){
            // Return iterator to lower/upper bound
            size_t low, up;
            low = std::distance( negScoreAll_.begin(),
                                 std::lower_bound( negScoreAll_.begin(),
                                                   negScoreAll_.end(),
                                                   posScoreAll_[i] ) );
            up = std::distance( negScoreAll_.begin(),
                                std::upper_bound( negScoreAll_.begin(),
                                                  negScoreAll_.end(),
                                                  posScoreAll_[i] ) );
            float p = 1.0f - ( float )( up + low ) / ( 2.0f * ( float ) negScoreAll_.size() );
            // avoid the rounding errors, such as p-value = 0 or p-value > 1
            if( p < 1.e-6 ) p = 1.e-6;
            if( p > 1.0f )  p = 1.0f;
            MOPS_Pvalue_.push_back( p );
        }
    }

    // for ZOOPS model:
    if( zoops_ ){
        // Sort log odds scores in ascending order
        std::sort( negScoreMax_.begin(), negScoreMax_.end(), std::less<float>() );
        std::sort( posScoreMax_.begin(), posScoreMax_.end(), std::less<float>() );
        for( size_t i = 0; i < posScoreMax_.size(); i++ ){
            // Return iterator to lower/upper bound
            size_t low, up;
            low = std::distance( negScoreMax_.begin(),
                                 std::lower_bound( negScoreMax_.begin(),
                                                   negScoreMax_.end(),
                                                   posScoreMax_[i] ) );
            up = std::distance( negScoreMax_.begin(),
                                std::upper_bound( negScoreMax_.begin(),
                                                  negScoreMax_.end(),
                                                  posScoreMax_[i] ) );
            float p = 1.0f - ( float )( up + low ) / ( 2.0f * ( float )negScoreMax_.size() );
            // avoid the rounding errors, such as p-value = 0 or p-value > 1
            if( p < 1.e-6 ) p = 1.e-6;
            if( p > 1.0f )  p = 1.0f;
            ZOOPS_Pvalue_.push_back( p );
        }
    }
}

void FDR::print(){

}

void FDR::write( char* odir, std::string basename ){

	std::string opath = std::string( odir ) + '/' + basename;

	if( savePRs_ ){
		/**
		 * save FDR results in two flat files for obtaining AUSFC:
		 * (1) posSequenceBasename.zoops.stats:
		 * TP, FP, FDR, recall, p-values, mFold and fractional occurrence for ZOOPS model
		 * (2) posSequenceBasename.mops.stats:
		 * TP, FP, FDR, recall and multiple occurrence for MOPS model
		 */

		// for ZOOPS model:
		if( zoops_ ){
			std::string opath_zoops_stats = opath + ".zoops.stats";
			std::ofstream ofile_zoops( opath_zoops_stats );
			// the headers:
			ofile_zoops << "TP" 	<< '\t'
						<< "FP" 	<< '\t'
						<< "FDR" 	<< '\t'
						<< "Recall"	<< '\t'
						<< "p-value"<< '\t'
						<< ( float )negSeqs_.size() / ( float )posSeqs_.size() << '\t'
                        << occ_frac_ << std::endl;

			for( size_t i = 0; i < ZOOPS_FDR_.size(); i++ ){
				ofile_zoops << ZOOPS_TP_[i]  << '\t'
							<< ZOOPS_FP_[i]  << '\t'
							<< ZOOPS_FDR_[i] << '\t'
							<< ZOOPS_Rec_[i] << '\t'
							<< PN_Pvalue_[i] << '\t'
							<< std::endl;
			}
		}

		// for MOPS model:
		if( mops_ ){
			std::string opath_mops_stats = opath + ".mops.stats";
			std::ofstream ofile_mops( opath_mops_stats );
			ofile_mops  << "TP" 	    << '\t'
						<< "FP" 	    << '\t'
						<< "FDR" 	    << '\t'
						<< "Recall"	    << '\t'
                        << occ_mult_    << std::endl;

			for( size_t i = 0; i < MOPS_FDR_.size(); i++ ){
				ofile_mops  << MOPS_TP_[i]  << '\t'
							<< MOPS_FP_[i]  << '\t'
							<< MOPS_FDR_[i] << '\t'
							<< MOPS_Rec_[i] << '\t' << std::endl;
			}
		}
	}

	if( savePvalues_ ){
		/**
		 * save FDR results in two flat files for ranking motifs:
		 * (1) posSequenceBasename.zoops.pvalues:	p-values for ZOOPS model
		 * (2) posSequenceBasename.mops.pvalues:	p-values for MOPS model
		 */
		if( zoops_ ){
			std::string opath_zoops = opath + ".zoops.pvalues";
			std::ofstream ofile_zoops( opath_zoops );
			for( size_t i = 0; i < ZOOPS_Pvalue_.size(); i++ ){
				ofile_zoops << std::setprecision( 3 ) << ZOOPS_Pvalue_[i] << std::endl;
			}
		}

		if( mops_ ){
			std::string opath_mops = opath + ".mops.pvalues";
			std::ofstream ofile_mops( opath_mops );
			for( size_t i = 0; i < MOPS_Pvalue_.size(); i++ ){
				ofile_mops  << std::setprecision( 3 ) << MOPS_Pvalue_[i] << std::endl;
			}
		}
	}

	if( saveLogOdds_ ){

		/**
		 * save log odds scores into two files for
		 * plotting the distribution of log odds scores:
		 * (1) posSequenceBasename.zoops.logOdds
		 * (2) posSequenceBasename.mops.logOdds
		 */
		if( zoops_ ){
			std::string opath_zoops_logOdds = opath + ".zoops.logOdds";
			std::ofstream ofile_zoops_logOdds( opath_zoops_logOdds );
			ofile_zoops_logOdds << "positive" << '\t'
								<< "negative" << std::endl;
			for( size_t i = 0; i < posScoreMax_.size(); i++ ){
				ofile_zoops_logOdds	<< std::setprecision( 6 )
									<< posScoreMax_[i] << '\t'
									<< negScoreMax_[i*negSeqs_.size()/posSeqs_.size()]
									<< std::endl;
			}
		}

		if( mops_ ){
			std::string opath_mops_logOdds = opath + ".mops.logOdds";
			std::ofstream ofile_mops_logOdds( opath_mops_logOdds );
			ofile_mops_logOdds 	<< "positive" << '\t'
								<< "negative" << std::endl;
			for( size_t i = 0; i < posScoreAll_.size(); i++ ){
				ofile_mops_logOdds 	<< std::setprecision( 6 )
									<< posScoreAll_[i] << '\t'
									<< negScoreAll_[i*negSeqs_.size()/posSeqs_.size()]
									<< std::endl;
			}
		}
	}
}
