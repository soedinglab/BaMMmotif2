#include <iomanip>
#include <chrono>

#include "EM.h"
#include "GibbsSampling.h"
#include "../seq_generator/SeqGenerator.h"

int main( int nargs, char* args[] ){

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

            if( Global::optimizePos and Global::savePi ) {
                // pre-define hyper-parameters for optimizing position priors
                size_t LW1 = posSet[0]->getL()-motif->getW()+1;
                std::string opath = std::string( Global::outputDirectory ) +'/'
                                    + Global::outputFileBasename + "_motif_"
                                    + std::to_string(n + 1) + ".pi";
                std::ofstream ofile( opath.c_str() );;
                for (size_t i = 1; i <= LW1; i++) {
                    ofile << model.getPi()[i] << std::endl;
                }
            }

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

        if (motif) delete motif;
    }

	return 0;
}
