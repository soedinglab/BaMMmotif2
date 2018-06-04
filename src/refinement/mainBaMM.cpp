#include <iomanip>
#include <chrono>

#include "Global.h"
#include "EM.h"
#include "GibbsSampling.h"
#include "../evaluation/FDR.h"

int main( int nargs, char* args[] ){

    auto t0_wall = std::chrono::high_resolution_clock::now();

    std::cout << std::endl
              << "======================================" << std::endl
              << "=      Welcome to use BaMM!motif     =" << std::endl
              << "=                   Version 2.0      =" << std::endl
              << "=              by Soeding Group      =" << std::endl
              << "=  http://www.mpibpc.mpg.de/soeding  =" << std::endl
              << "======================================" << std::endl;

	// seed random number
	srand( 42 );
	Global::rngx.seed( 42 );

	// initialization
	Global::init( nargs, args );

	if( Global::verbose ){
		std::cout << std::endl
                  << "************************" << std::endl
                  << "*   Background Model   *" << std::endl
                  << "************************" << std::endl;
	}

    std::vector<Sequence*> posSet = Global::posSequenceSet->getSequences();

	BackgroundModel* bgModel;
	if( !Global::bgModelGiven ){
		bgModel = new BackgroundModel( posSet,
                                       Global::bgModelOrder,
                                       Global::bgModelAlpha,
                                       Global::interpolateBG,
                                       Global::outputFileBasename );

    } else {
		bgModel = new BackgroundModel( Global::bgModelFilename );

	}

    // always save the background model
    bgModel->write( Global::outputDirectory, Global::outputFileBasename );

	if( Global::verbose ){
        std::cout << std::endl
                  << "***************************" << std::endl
                  << "*   Initial Motif Model   *" << std::endl
                  << "***************************" << std::endl;
	}

	MotifSet motif_set( Global::initialModelFilename,
                        Global::addColumns.at(0),
                        Global::addColumns.at(1),
                        Global::initialModelTag,
                        Global::posSequenceSet,
                        Global::posSequenceSet->getBaseFrequencies(),
                        Global::modelOrder,
                        Global::modelAlpha,
                        Global::maxPWM,
                        Global::q );

    /**
     * Filter out short sequences
     */
    std::vector<Sequence*>::iterator it = posSet.begin();
    while( it != posSet.end() ){
        if( (*it)->getL() < motif_set.getMaxW() ){
            //std::cout << "Warning: remove the short sequence: " << (*it)->getHeader() << std::endl;
            posSet.erase(it);
        } else {
            it++;
        }
    }

	if( Global::verbose ){
        std::cout << std::endl
                  << "*********************" << std::endl
                  << "*   BaMM Training   *" << std::endl
                  << "*********************" << std::endl;
	}

    // sample negative sequence set B1set based on s-mer frequencies
    // from positive training sequence set
    std::vector<Sequence*>  negSet;
    std::vector<std::unique_ptr<Sequence>> negSeqs;
    size_t minSeqN = 5000;
    bool rest = minSeqN % posSet.size();
    if( posSet.size() < minSeqN ){
        Global::mFold = minSeqN / posSet.size() + rest;
    }

    SeqGenerator negseq( posSet, NULL, Global::sOrder, 1.0f, Global::genericNeg );
    negSeqs = negseq.sample_bgseqset_by_fold( Global::mFold );

    // convert unique_ptr to regular pointer
    for( size_t n = 0; n < negSeqs.size(); n++ ) {
//        negSeqs[n]->print();    // print out negative sequences for checking
        negSet.push_back( negSeqs[n].release() );
        negSeqs[n].get_deleter();
    }

//#pragma omp parallel for
    for( size_t n = 0; n < motif_set.getN(); n++ ){
		// deep copy each motif in the motif set
		Motif* motif = new Motif( *motif_set.getMotifs()[n] );

		if( Global::saveInitialBaMMs ){
			// optional: save initial model
			motif->write( Global::outputDirectory,
                          Global::outputFileBasename + "_init_motif_" + std::to_string( n+1 ) );
		}

		// optimize the model with either EM or Gibbs sampling
		if( Global::EM ){
			EM model( motif, bgModel, posSet, Global::optimizeQ, Global::verbose, Global::f );
			// learn motifs by EM
			if( !Global::advanceEM ) {
                model.optimize();
            } else {
                model.mask();
            }

            // write model parameters on the disc
			if( Global::saveBaMMs ){
				model.write( Global::outputDirectory,
                             Global::outputFileBasename + "_motif_" + std::to_string( n+1 ),
                             Global::ss );
			}

            // print out the optimized q for checking:
            std::cout << "optimized q = " << model.getQ() << std::endl;

		} else if ( Global::CGS ){
			GibbsSampling model( motif, bgModel, posSet,
                                 !Global::noQSampling, Global::verbose );
			// learn motifs by collapsed Gibbs sampling
			model.optimize();
			// write model parameters on the disc
			if( Global::saveBaMMs ){
				model.write( Global::outputDirectory,
                             Global::outputFileBasename + "_motif_" + std::to_string( n+1 ),
                             Global::ss );
			}

            // print out the optimized q for checking:
            std::cout << "optimized q = " << model.getQ() << std::endl;

		} else {
			std::cout << "Note: the model is not optimized!\n";
		}

        // write out the learned model
        motif->write( Global::outputDirectory,
                      Global::outputFileBasename + "_motif_" + std::to_string( n+1 ) );

        if( Global::scoreSeqset ){
            // score the model on sequence set
            if( Global::verbose ) {
                std::cout << std::endl
                          << "*************************" << std::endl
                          << "*    Score Sequences    *" << std::endl
                          << "*************************" << std::endl
                          << std::endl;
            }

            // Define bg model depending on motif input, and learning
            // use bgModel generated from input sequences when prediction is turned on
            BackgroundModel* bg = bgModel;
            // if no optimization is applied, get bgModel from the input
            if( !Global::EM and !Global::CGS ) {
                // use provided bgModelFile if initialized with bamm format
                if( Global::initialModelTag == "BaMM" ) {
                    if( Global::bgModelFilename == NULL ) {
                        std::cout << "No background Model file provided for initial search motif!\n";
                        exit( 1 );
                    }
                    bg = new BackgroundModel( Global::bgModelFilename );
                } else if( Global::initialModelTag == "PWM" ){
                    // this means that also the global motif order needs to be adjusted;
                    Global::modelOrder = 0;
                }
            }

            ScoreSeqSet scoreNegSet( motif, bgModel, negSet );
            scoreNegSet.calcLogOdds();

            // print out log odds scores for checking before reranking
            if( Global::saveLogOdds ){
                scoreNegSet.writeLogOdds(Global::outputDirectory,
                                         Global::outputFileBasename + ".negSet",
                                         Global::ss );
            }

            std::vector<std::vector<float>> negAllScores = scoreNegSet.getMopsScores();
            std::vector<float> negScores;
            for( size_t n = 0; n < negSet.size(); n++ ){
                negScores.insert( std::end( negScores ),
                                  std::begin( negAllScores[n] ),
                                  std::end( negAllScores[n] ) );
            }

            // calculate p-values based on positive and negative scores
            ScoreSeqSet scorePosSet( motif, bg, posSet );
            scorePosSet.calcLogOdds();

            // print out log odds scores for checking before reranking
            if( Global::saveLogOdds ){
                scorePosSet.writeLogOdds(Global::outputDirectory,
                                         Global::outputFileBasename + "_motif_" + std::to_string( n+1 ),
                                         Global::ss );
            }

            std::vector<std::vector<float>> posScores = scorePosSet.getMopsScores();
            scorePosSet.calcPvalues( posScores, negScores );

            scorePosSet.write( Global::outputDirectory,
                               Global::outputFileBasename + "_motif_" + std::to_string( n+1 ),
                               Global::pvalCutoff,
                               Global::ss );

        }

        if( motif )		delete motif;
	}

    // evaluate motifs
	if( Global::FDR ){
		if( Global::verbose ){
            std::cout << std::endl
                      << "***********************" << std::endl
                      << "*   BaMM validation   *" << std::endl
                      << "***********************" << std::endl;
		}

        /**
         * cross-validate the motif model
         */
//#pragma omp parallel for
        for( size_t n = 0; n < motif_set.getN(); n++ ){
			Motif* motif = new Motif( *motif_set.getMotifs()[n] );
			FDR fdr( posSet, negSet,
                     motif, bgModel, Global::cvFold,
                     Global::mops, Global::zoops,
                     Global::savePRs, Global::savePvalues, Global::saveLogOdds );
			fdr.evaluateMotif( Global::EM, Global::CGS, Global::optimizeQ, Global::advanceEM, Global::f );
			fdr.write( Global::outputDirectory,
                       Global::outputFileBasename + "_motif_" + std::to_string( n+1 ) );
			if( motif )		delete motif;
		}

    }

    std::cout << std::endl
              << "******************" << std::endl
              << "*   Statistics   *" << std::endl
              << "******************" << std::endl;
	Global::printStat();

    auto t1_wall = std::chrono::high_resolution_clock::now();
    auto t_diff = std::chrono::duration_cast<std::chrono::duration<double>>(t1_wall-t0_wall);
    std::cout << std::endl << "------ Runtime: " << t_diff.count() <<" seconds -------" << std::endl;

	// free memory
	if( bgModel ) delete bgModel;

	Global::destruct();

	return 0;
}
