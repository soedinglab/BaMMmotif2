//
// Created by wanwan on 04.09.17.
// This function is used for scanning motif occurrences in a sequence set
//

#include "ScoreSeqSet.h"
#include "../init/MotifSet.h"
#include "../seq_generator/SeqGenerator.h"

int main( int nargs, char* args[] ) {

    /**
     * initialization
     */
	Global G( nargs, args );

    /**
     * Build up the background model
     */
    // use bgModel generated from input sequences when prediction is turned on
    BackgroundModel* bgModel = new BackgroundModel( Global::negSequenceSet->getSequences(),
                                                    Global::outputFileBasename );

    // use provided bgModelFile if initialized with bamm format
    if( Global::initialModelTag == "BaMM" ) {
        if( Global::bgModelFilename == NULL ) {
            std::cerr << "Error: No background model file provided for initial "
                    "search motif!" << std::endl;
            exit( 1 );
        }
        // get background model from the given file
        bgModel = new BackgroundModel( Global::bgModelFilename );
    } else if( Global::initialModelTag == "PWM" ){
        // this means that also the global motif order needs to be adjusted;
        Global::modelOrder = 0;
    }

    if( Global::saveBgModel ){
        // save background model
        bgModel->write( Global::outputDirectory, Global::outputFileBasename );
    }

    /**
     * Initialize the model
     */
    MotifSet motif_set( Global::initialModelFilename,
                        Global::addColumns.at(0),
                        Global::addColumns.at(1),
                        Global::initialModelTag,
                        Global::posSequenceSet );

    std::vector<Sequence*> posSet = Global::posSequenceSet->getSequences();

    /**
     * Filter out short sequences
     */
    auto it = posSet.begin();
    while( it != posSet.end() ){
        if( (*it)->getL() < motif_set.getMaxW() ){
            std::cout << "Warning: remove the short sequence: "
                      << (*it)->getHeader()
                      << std::endl;
            posSet.erase(it);
        } else {
            it++;
        }
    }

    std::vector<Sequence *> negSet;
    if( !Global::noNegset ) {
        /**
         * Sample negative sequence set based on s-mer frequencies
         */
        // sample negative sequence set B1set based on s-mer frequencies
        // from positive training sequence set
        std::vector<std::unique_ptr<Sequence>> negSeqs;
        SeqGenerator negseq(posSet);
        size_t posN = posSet.size();

        size_t mFold = Global::mFold;
        if (posN * mFold < Global::minNegN) {
            auto rest = bool(Global::minNegN % posN);
            mFold = Global::minNegN / posN + rest;
        }

        negSeqs = negseq.sample_bgset_by_fold(mFold);
        if (Global::verbose) {
            std::cout << "\n" << negSeqs.size()
                      << " background sequences are generated."
                      << std::endl;
        }

        // convert unique_ptr to regular pointer
        size_t negN = negSeqs.size();
        for (size_t n = 0; n < negN; n++) {
            negSet.push_back(negSeqs[n].release());
            negSeqs[n].get_deleter();
        }
    }

//#pragma omp parallel for
    for( size_t n = 0; n < motif_set.getN(); n++ ) {
        // deep copy each motif in the motif set
        Motif *motif = new Motif( *motif_set.getMotifs()[n] );

        std::string fileExtension;
        if( Global::initialModelTag == "PWM" ){
            fileExtension = "_motif_" + std::to_string( n+1 );
        }

        if( Global::saveInitialBaMMs ){
            // write out the foreground model
            motif->write( Global::outputDirectory,
                          Global::outputFileBasename + fileExtension );
        }

        // score positive sequence set
        // calculate p-values based on positive and negative scores
        ScoreSeqSet scorePosSet( motif, bgModel, posSet );
        scorePosSet.calcLogOdds();


        if( !Global::noNegset ) {
            // score negative sequence set
            ScoreSeqSet scoreNegSet(motif, bgModel, negSet);
            scoreNegSet.calcLogOdds();
            std::vector<std::vector<float>> negAllScores = scoreNegSet.getMopsScores();
            std::vector<float> negScores;
            for (size_t i = 0; i < negSet.size(); i++) {
                negScores.insert(std::end(negScores),
                                 std::begin(negAllScores[i]),
                                 std::end(negAllScores[i]));
            }
            std::vector<std::vector<float>> posScores = scorePosSet.getMopsScores();
            scorePosSet.calcPvalues( posScores, negScores );
        }

        scorePosSet.write( Global::outputDirectory,
                           Global::outputFileBasename + fileExtension,
                           Global::pvalCutoff,
                           Global::logOddsCutoff,
                           Global::ss );

        delete motif;
    }

    delete bgModel;

    return 0;
}
