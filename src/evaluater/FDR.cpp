#include "FDR.h"

FDR::FDR( std::vector<Sequence*> posSeqs, std::vector<Sequence*> negSeqs,
          Motif* motif, BackgroundModel* bgModel){

	posSeqs_	= posSeqs;
	negSeqs_	= negSeqs;

	motif_ 		= motif;
    bgModel_    = bgModel;
    q_          = motif_->getQ();

	occ_frac_	= 0.0f;
	occ_mult_	= 0.0f;

}

FDR::~FDR(){

}

void FDR::evaluateMotif( size_t perLoopThreads ){

	std::vector<std::vector<float>> mops_scores;
	std::vector<float> 				zoops_scores;
    float updatedQ = q_;            // obtain the updated q

    /**
	 * Cross validation
	 */
#pragma omp parallel for num_threads(perLoopThreads)
	for( size_t fold = 0; fold < Global::cvFold; fold++ ){

		// deep copy the initial motif
		Motif* motif = new Motif( *motif_ );

		/**
		 * Draw sequences for each training, test and negative sets
		 */
		std::vector<Sequence*> testSet;
		std::vector<Sequence*> trainSet;
        std::vector<Sequence*> negSet;
		for( size_t n = 0; n <= posSeqs_.size() - 1; n += Global::cvFold ){
			for( size_t f = 0; f < Global::cvFold; f++ ){
				if( n+f < posSeqs_.size() ){
                    if( f != fold ){
                        trainSet.push_back( posSeqs_[n+f] );
                    } else {
                        testSet.push_back( posSeqs_[n+f] );
                    }
                }
			}
		}
		for( size_t n = 0; n <= negSeqs_.size()-Global::cvFold; n += Global::cvFold ){
            negSet.push_back( negSeqs_[n] );
		}

		/**
		 * Training
		 */
		// learn motif from each training set
		if( Global::EM ){
			EM model( motif, bgModel_, trainSet );
            if( Global::advanceEM ){
                model.mask();
            } else {
                model.optimize();
            }
            updatedQ = model.getQ();
		} else if ( Global::CGS ){
			GibbsSampling model( motif, bgModel_, trainSet );
			model.optimize();
            updatedQ = model.getQ();
		}

		/**
		 * Testing
		 */
		// score positive test sequences with (learned) motif
		ScoreSeqSet score_testset( motif, bgModel_, testSet );
        score_testset.calcLogOdds();

        // score negative sequence set
        ScoreSeqSet score_negset( motif, bgModel_, negSet );
        score_negset.calcLogOdds();

#pragma omp critical
        {
            if (Global::mops) {
                mops_scores = score_testset.getMopsScores();

                for (size_t n = 0; n < testSet.size(); n++) {
                    posScoreAll_.insert(std::end(posScoreAll_),
                                        std::begin(mops_scores[n]),
                                        std::end(mops_scores[n]));
                }

                mops_scores.clear(); // reuse mops_scores vector

                mops_scores = score_negset.getMopsScores();
                for (size_t n = 0; n < negSet.size(); n++) {
                    negScoreAll_.insert(std::end(negScoreAll_),
                                        std::begin(mops_scores[n]),
                                        std::end(mops_scores[n]));
                }
            }

            if (Global::zoops) {
                zoops_scores = score_testset.getZoopsScores();
                posScoreMax_.insert(std::end(posScoreMax_),
                                    std::begin(zoops_scores),
                                    std::end(zoops_scores));

                zoops_scores.clear(); // reuse zoops_scores vector
                zoops_scores = score_negset.getZoopsScores();
                negScoreMax_.insert(std::end(negScoreMax_),
                                    std::begin(zoops_scores),
                                    std::end(zoops_scores));
            }

        }
		if( motif ) 				delete motif;
	}

    // update Q
    q_ = updatedQ;

    // calculate precision and recall
    calculatePR();

    if( Global::savePvalues ){

        std::cout << " ____________________" << std::endl;
        std::cout << "|                    |" << std::endl;
        std::cout << "| calculate P-values |" << std::endl;
        std::cout << "|____________________|" << std::endl;
        std::cout << std::endl;

		calculatePvalues();
	}
}

void FDR::calculatePR(){

	size_t posN = posSeqs_.size();
	size_t negN = negSeqs_.size();

	float mFold = ( float )negN / ( float )posN;

    srand(42);

	// for MOPS model:
	if( Global::mops ){
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
	if( Global::zoops ){
		// Sort log odds scores in descending order
		std::sort( posScoreMax_.begin(), posScoreMax_.end(), std::greater<float>() );
		std::sort( negScoreMax_.begin(), negScoreMax_.end(), std::greater<float>() );

		// Rank and score these log odds score values
		size_t idx_posMax = 0;
		size_t idx_negMax = 0;
		size_t min_idx_pos = 0;
		size_t posN_est = size_t( q_ * ( float )posN );

        // set limit for using the exponential extrapolation for p-value calculation
        size_t n_top = ( size_t )std::fmin(100, negN / 10);

        float lambda = Global::epsilon;
        for( size_t l = 0; l < n_top; l++ ){
            lambda += negScoreMax_[l] - negScoreMax_[n_top];
        }
        lambda /= n_top;

        float Sl;

        // dynamic programming for scoring both positive and negative sets
		for( size_t i = 0; i < posN + negN; i++ ){

            if( idx_posMax == 0 || idx_negMax == negN ){
                Sl = posScoreMax_[idx_posMax];
                idx_posMax++;
            } else if( idx_posMax == posN ){
                Sl = negScoreMax_[idx_negMax];
                idx_negMax++;
            } else if( idx_posMax < posN and idx_negMax < negN ) {
                if (posScoreMax_[idx_posMax] > negScoreMax_[idx_negMax]) {
                    Sl = posScoreMax_[idx_posMax];
                    idx_posMax++;
                } else if (posScoreMax_[idx_posMax] == negScoreMax_[idx_negMax]
                           and rand() % 2 == 0) {
                    Sl = posScoreMax_[idx_posMax];
                    idx_posMax++;
                } else {
                    Sl = negScoreMax_[idx_negMax];
                    idx_negMax++;
                }
            }

            // calculate p-values in two different ways:
            float p_value;
            if( Sl <= negScoreMax_[n_top] ){
                float Sl_upper = *(std::lower_bound( negScoreMax_.begin(),
                                                     negScoreMax_.end(), Sl,
                                                     std::greater<float>() )-1);
                float Sl_lower = *std::upper_bound( negScoreMax_.begin(),
                                                    negScoreMax_.end(), Sl,
                                                    std::greater<float>() );
                p_value = (idx_negMax + ( Sl_upper- Sl)
                                        / (Sl_upper - Sl_lower + Global::epsilon))
                          / (float)negN;

            } else {
                // p-value is calculated by relying on a parametric fit
                // of the exponentially cumulative distribution
                p_value = n_top * expf( ( negScoreMax_[n_top] - Sl ) / lambda ) / negN;
            }
            if( p_value < Global::epsilon ) p_value = Global::epsilon;
            if( p_value > 1.0f )            p_value = 1.0f;
            //assert( p_value >= 0.f && p_value <= 1.f );

			PN_Pvalue_.push_back( p_value );

            // calculate TP, FP, FDR and Recall
            float TP = ( float )idx_posMax;
            float FP = ( float )idx_negMax / mFold;
            ZOOPS_TP_.push_back( TP );
            ZOOPS_FP_.push_back( FP );
            ZOOPS_FDR_.push_back( FP / ( TP + FP ) );
            ZOOPS_Rec_.push_back( TP / ( float )posN );

			// take the faction of q sequences as real positives
			if( idx_posMax == posN_est ){
                min_idx_pos = i;
			}
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
    if( Global::mops ){
        // Sort log odds scores in ascending order
        std::sort( negScoreAll_.begin(), negScoreAll_.end(), std::less<float>() );
        std::sort( posScoreAll_.begin(), posScoreAll_.end(), std::less<float>() );
        for( size_t i = 0; i < posScoreAll_.size(); i++ ){
            // Return iterator to lower/upper bound
            size_t low = ( size_t )std::distance( negScoreAll_.begin(),
                                                  std::lower_bound( negScoreAll_.begin(),
                                                                    negScoreAll_.end(),
                                                                    posScoreAll_[i] ) );
            size_t up = ( size_t )std::distance( negScoreAll_.begin(),
                                                 std::upper_bound( negScoreAll_.begin(),
                                                                   negScoreAll_.end(),
                                                                   posScoreAll_[i] ) );
            float p = 1.0f - ( float )( up + low ) / ( 2.0f * ( float ) negScoreAll_.size() );
            // avoid the rounding errors, such as p-value = 0 or p-value > 1
            if( p < Global::epsilon )   p = Global::epsilon;
            if( p > 1.0f )              p = 1.0f;
            MOPS_Pvalue_.push_back( p );
        }
    }

    // for ZOOPS model:
    if( Global::zoops ){
        // Sort log odds scores in ascending order
        std::sort( negScoreMax_.begin(), negScoreMax_.end(), std::less<float>() );
        std::sort( posScoreMax_.begin(), posScoreMax_.end(), std::less<float>() );
        for( size_t i = 0; i < posScoreMax_.size(); i++ ){
            // Return iterator to lower/upper bound
            size_t low = ( size_t )std::distance( negScoreMax_.begin(),
                                                  std::lower_bound( negScoreMax_.begin(),
                                                                    negScoreMax_.end(),
                                                                    posScoreMax_[i] ) );
            size_t up = ( size_t )std::distance( negScoreMax_.begin(),
                                                 std::upper_bound( negScoreMax_.begin(),
                                                                   negScoreMax_.end(),
                                                                   posScoreMax_[i] ) );
            float p = 1.0f - ( float )( up + low ) / ( 2.0f * ( float )negScoreMax_.size() );
            // avoid the rounding errors, such as p-value = 0 or p-value > 1
            if( p < Global::epsilon )   p = Global::epsilon;
            if( p > 1.0f )              p = 1.0f;
            ZOOPS_Pvalue_.push_back( p );
        }
    }
}

void FDR::print(){

}

void FDR::write( char* odir, std::string basename ){

	std::string opath = std::string( odir ) + '/' + basename;

	if( Global::savePRs ){
		/**
		 * save FDR results in two flat files for obtaining AvRec:
		 * (1) posSequenceBasename.zoops.stats:
		 * TP, FP, FDR, recall, p-values, mFold and fractional occurrence for ZOOPS model
		 * (2) posSequenceBasename.mops.stats:
		 * TP, FP, FDR, recall and multiple occurrence for MOPS model
		 */

		// for ZOOPS model:
		if( Global::zoops ){
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
		if( Global::mops ){
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

	if( Global::savePvalues ){
		/**
		 * save FDR results in two flat files for ranking motifs:
		 * (1) posSequenceBasename.zoops.pvalues:	p-values for ZOOPS model
		 * (2) posSequenceBasename.mops.pvalues:	p-values for MOPS model
		 */
		if( Global::zoops ){
			std::string opath_zoops = opath + ".zoops.pvalues";
			std::ofstream ofile_zoops( opath_zoops );
			for( size_t i = 0; i < ZOOPS_Pvalue_.size(); i++ ){
				ofile_zoops << std::setprecision( 3 ) << ZOOPS_Pvalue_[i] << std::endl;
			}
		}

		if( Global::mops ){
			std::string opath_mops = opath + ".mops.pvalues";
			std::ofstream ofile_mops( opath_mops );
			for( size_t i = 0; i < MOPS_Pvalue_.size(); i++ ){
				ofile_mops  << std::setprecision( 3 ) << MOPS_Pvalue_[i] << std::endl;
			}
		}
	}

	if( Global::saveLogOdds ){
		/**
		 * save log odds scores into two files for
		 * plotting the distribution of log odds scores:
		 * (1) posSequenceBasename.zoops.logOdds
		 * (2) posSequenceBasename.mops.logOdds
		 */

		if( Global::zoops ){
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

		if( Global::mops ){
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

void FDR::saveUnsortedLogOdds( std::string opath, std::vector<float> logOdds ){

    std::ofstream ofile( opath );

    for( size_t i = 0; i < logOdds.size(); i++ ){
        ofile << i+1 << '\t' << std::setprecision( 6 ) << logOdds[i] << std::endl;
    }
}
