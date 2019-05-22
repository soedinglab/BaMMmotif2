#include <iomanip>
#include <chrono>

#include "../Global/Global.h"
#include "EM.h"
#include "GibbsSampling.h"
#include "../evaluater/FDR.h"

int main( int nargs, char* args[] ){

    auto t0_wall = std::chrono::high_resolution_clock::now();
    std::cout.precision(4);

    std::cout << std::endl
              << "======================================" << std::endl
              << "=      Welcome to use BaMM!motif     =" << std::endl
              << "=                    Version 2.0     =" << std::endl
              << "=               by Soeding Group     =" << std::endl
              << "=  http://www.mpibpc.mpg.de/soeding  =" << std::endl
              << "======================================" << std::endl;

    /**
     * initialization
     */
    Global G( nargs, args );

    /**
     * Build up the background model
     */
	BackgroundModel* bgModel;
	if( !Global::bgModelGiven ){
		bgModel = new BackgroundModel( Global::posSequenceSet->getSequences(),
                                       Global::outputFileBasename );
    } else {
		bgModel = new BackgroundModel( Global::bgModelFilename );
	}

    // always save the background model
    bgModel->write( Global::outputDirectory, Global::outputFileBasename );

    if( Global::verbose ){
        bgModel->print();
    }

    /**
     * Initialize the model
     */
	MotifSet motif_set( Global::initialModelFilename,
                        Global::addColumns.at(0),
                        Global::addColumns.at(1),
                        Global::initialModelTag,
                        Global::posSequenceSet );

    if( Global::verbose ){
        motif_set.print();
    }

    /**
     * Filter out short sequences
     */
    std::vector<Sequence*> posSet = Global::posSequenceSet->getSequences();
    std::vector<Sequence*>::iterator it = posSet.begin();
    while( it != posSet.end() ){
        if( (*it)->getL() < motif_set.getMaxW() ){
            if( Global::verbose ){
                std::cout << "Warning: remove the short sequence: "
                          << (*it)->getHeader() << std::endl;
            }
            posSet.erase(it);
        } else {
            it++;
        }
    }

    /**
     * Generate negative sequence set
     * for training, scoring and cross-validation
     */
    std::vector<Sequence*> negSet;
    if( Global::negSeqGiven ){
        // take the given negative sequence set
        negSet = Global::negSequenceSet->getSequences();
    } else {
        // generate negative sequence set based on s-mer frequencies
        // from positive training sequence set
        std::vector<std::unique_ptr<Sequence>> negSeqs;
        SeqGenerator negseq(posSet, nullptr);
        negSeqs = negseq.sample_bgset_by_fold(Global::mFold);

        // convert unique_ptr to regular pointer
        for (size_t n = 0; n < negSeqs.size(); n++) {
            negSet.push_back(negSeqs[n].release());
            negSeqs[n].get_deleter();
        }

    }

    /**
     * Train the model(s)
     */
    for( size_t n = 0; n < motif_set.getN(); n++ ) {
        // deep copy each motif in the motif set
        Motif *motif = new Motif(*motif_set.getMotifs()[n]);

        if (Global::saveInitialBaMMs) {
            // optional: save initial model
            motif->write(Global::outputDirectory,
                         Global::outputFileBasename + "_init_motif_" + std::to_string(n + 1));
        }

        // optimize the model with either EM or Gibbs sampling
        if (Global::EM) {
            if (Global::verbose) {
                std::cout << " _________________" << std::endl
                          << "|*               *|" << std::endl
                          << "|   EM training   |" << std::endl
                          << "|*_______________*|" << std::endl
                          << std::endl;
            }
            EM model(motif, bgModel, posSet);

            if( Global::optimizeQ ) {
                std::cout << std::endl << "Before, fraction parameter q="
                          << std::setprecision(4) << model.getQ() << std::endl;
            }
            // learn motifs by EM
            if (!Global::advanceEM) {
                model.optimize();

/*                // todo: print out for checking
                if(Global::optimizePos){
                    std::string opath = std::string( Global::outputDirectory ) + '/'
                                        + Global::outputFileBasename + "_motif_" +
                                        std::to_string(n + 1)+".pi";
                    std::ofstream ofile( opath.c_str() );
                    size_t LW1 = posSet[0]->getL()-motif->getW() +1;
                    for(size_t i = 1; i <= LW1; i++ ){
                        ofile << model.getPi()[i] << std::endl;
                    }
                }*/

            } else {
                model.mask();
            }

            // write model parameters on the disc
            if (Global::saveBaMMs) {
                model.write(Global::outputDirectory,
                            Global::outputFileBasename + "_motif_" + std::to_string(n + 1));
            }

            if (Global::verbose) {
                std::cout << " ____________________" << std::endl
                          << "|*                  *|" << std::endl
                          << "|   After training   |" << std::endl
                          << "|*__________________*|" << std::endl
                          << std::endl;

                model.print();
            }
            std::cout << std::endl << "After, fraction parameter q="
                      << std::setprecision(4) << model.getQ() << std::endl;

        } else if (Global::CGS) {

            if (Global::verbose) {
                std::cout << " ____________________" << std::endl
                          << "|*                  *|" << std::endl
                          << "|   Gibbs sampling   |" << std::endl
                          << "|*__________________*|" << std::endl
                          << std::endl;
            }

            GibbsSampling model(motif, bgModel, posSet);

            // learn motifs by collapsed Gibbs sampling
            model.optimize();
            // write model parameters on the disc
            if (Global::saveBaMMs) {
                model.write(Global::outputDirectory,
                            Global::outputFileBasename + "_motif_" + std::to_string(n + 1),
                            Global::ss);
            }

            if( Global::optimizeQ ) {
                // print out the optimized q for checking:
                std::cout << "optimized q = " << model.getQ() << std::endl;
            }
        } else {
            std::cout << "\nNote: the model is not optimized!\n";
        }

        // write out the (learned) foreground model
        motif->write(Global::outputDirectory,
                     Global::outputFileBasename + "_motif_" + std::to_string(n + 1));

        if (Global::scoreSeqset) {
            /**
             *  score the model on sequence set
             */

            // Define bg model depending on motif input, and learning
            // use bgModel generated from input sequences when prediction is turned on
            BackgroundModel *bg = bgModel;
            // if no optimization is applied, get bgModel from the input
            if (!Global::EM and !Global::CGS) {
                // use provided bgModelFile if initialized with bamm format
                if (Global::initialModelTag == "BaMM") {
                    if (Global::bgModelFilename == NULL) {
                        std::cout << "No background Model file provided for initial search motif!\n";
                        exit(1);
                    }
                    bg = new BackgroundModel(Global::bgModelFilename);
                } else if (Global::initialModelTag == "PWM") {
                    // this means that also the global motif order needs to be adjusted;
                    Global::modelOrder = 0;
                }
            }

            ScoreSeqSet scoreNegSet(motif, bg, negSet);

            scoreNegSet.calcLogOdds();

            // print out log odds scores for checking before re-ranking
            if (Global::saveLogOdds) {
                scoreNegSet.writeLogOdds(Global::outputDirectory,
                                         Global::outputFileBasename + ".negSet",
                                         Global::ss);
            }

            std::vector<std::vector<float>> negAllScores = scoreNegSet.getMopsScores();
            std::vector<float> negScores;
            for (size_t n = 0; n < negSet.size(); n++) {
                negScores.insert(std::end(negScores),
                                 std::begin(negAllScores[n]),
                                 std::end(negAllScores[n]));
            }

            // calculate p-values based on positive and negative scores
            ScoreSeqSet scorePosSet(motif, bg, posSet);
            scorePosSet.calcLogOdds();

            // print out log odds scores for checking before re-ranking
            if (Global::saveLogOdds) {
                scorePosSet.writeLogOdds(Global::outputDirectory,
                                         Global::outputFileBasename + "_motif_" + std::to_string(n + 1),
                                         Global::ss);
            }

            std::vector<std::vector<float>> posScores = scorePosSet.getMopsScores();
            scorePosSet.calcPvalues(posScores, negScores);

            // save the occurrences that has a p-value above certain cutoff
            scorePosSet.write(Global::outputDirectory,
                              Global::outputFileBasename + "_motif_" + std::to_string(n + 1),
                              Global::pvalCutoff,
                              Global::ss);

        }

        if (motif) delete motif;
    }

    // evaluate motifs
	if( Global::FDR ){
        /**
         * cross-validate the motif model
         */
//#pragma omp parallel for
        for( size_t n = 0; n < motif_set.getN(); n++ ){
			Motif* motif = new Motif( *motif_set.getMotifs()[n] );
			FDR fdr( posSet, negSet, motif, bgModel );
			fdr.evaluateMotif();
			fdr.write( Global::outputDirectory,
                       Global::outputFileBasename + "_motif_" + std::to_string( n+1 ) );
			if( motif )		delete motif;
		}

    }

    auto t1_wall = std::chrono::high_resolution_clock::now();
    auto t_diff = std::chrono::duration_cast<std::chrono::duration<double>>(t1_wall-t0_wall);
    std::cout << std::endl << "---- Total Runtime: " << t_diff.count() <<" seconds ----" << std::endl;

	// free memory
	if( bgModel ) delete bgModel;
    if( !Global::negSeqGiven and negSet.size() ){
        for (size_t n = 0; n < negSet.size(); n++) {
            delete negSet[n];
        }
    }

	return 0;
}
