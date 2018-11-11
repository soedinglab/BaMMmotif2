//
// Created by wanwan on 04.09.17.
// This function is used for scanning motif occurrences in a sequence set
//

#include "GScan.h"
#include "ScoreSeqSet.h"
#include "../init/MotifSet.h"
#include "../seq_generator/SeqGenerator.h"

int main( int nargs, char* args[] ) {

    /**
     * initialization
     */

	GScan::init( nargs, args );

    /**
     * Build up the background model
     */
    // use bgModel generated from input sequences when prediction is turned on
    BackgroundModel* bgModel = new BackgroundModel( GScan::negSequenceSet->getSequences(),
                                       GScan::bgModelOrder,
                                       GScan::bgModelAlpha,
                                       GScan::interpolateBG,
                                       GScan::outputFileBasename );

    // use provided bgModelFile if initialized with bamm format
    if( GScan::initialModelTag == "BaMM" ) {
        if( GScan::bgModelFilename == NULL ) {
            std::cerr << "Error: No background model file provided for initial search motif!" << std::endl;
            exit( 1 );
        }
        // get background model from the given file
        bgModel = new BackgroundModel( GScan::bgModelFilename );
    } else if( GScan::initialModelTag == "PWM" ){
        // this means that also the global motif order needs to be adjusted;
        GScan::modelOrder = 0;
    }

    if(GScan::saveInitialModel){
        // save background model
        bgModel->write(GScan::outputDirectory, GScan::outputFileBasename);
    }

    /**
     * Initialize the model
     */
    MotifSet motif_set( GScan::initialModelFilename,
                        GScan::addColumns.at(0),
                        GScan::addColumns.at(1),
                        GScan::initialModelTag,
                        GScan::posSequenceSet,
                        bgModel->getV(),
                        GScan::bgModelOrder,
                        GScan::modelOrder,
                        GScan::modelAlpha,
                        GScan::maxPWM);

    std::vector<Sequence*> posSet = GScan::posSequenceSet->getSequences();
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

    /**
     * Sample negative sequence set based on s-mer frequencies
     */
    std::vector<Sequence*>  negSet;
    size_t minSeqN = 5000;
    // sample negative sequence set B1set based on s-mer frequencies
    // from positive training sequence set
    std::vector<std::unique_ptr<Sequence>> negSeqs;
    SeqGenerator negseq( posSet );
    size_t posN = posSet.size();
    bool rest = minSeqN % posN;
    if ( posN * GScan::mFold < minSeqN ) {
        GScan::mFold = minSeqN / posN + rest;
    }
    negSeqs = negseq.sample_bgseqset_by_fold(GScan::mFold);
    std::cout << negSeqs.size() << " background sequences are generated." << std::endl;

    // convert unique_ptr to regular pointer
    for( size_t n = 0; n < negSeqs.size(); n++ ) {
        negSet.push_back( negSeqs[n].release() );
        negSeqs[n].get_deleter();
    }

#pragma omp parallel for
    for( size_t n = 0; n < motif_set.getN(); n++ ) {
        // deep copy each motif in the motif set
        Motif *motif = new Motif( *motif_set.getMotifs()[n] );

        std::string fileExtension;
        if( GScan::initialModelTag == "PWM" ){
            fileExtension = "_motif_" + std::to_string( n+1 );
        }

        if( GScan::saveInitialModel ){
            // write out the foreground model
            motif->write( GScan::outputDirectory,
                          GScan::outputFileBasename + fileExtension );
        }

        // score negative sequence set
        ScoreSeqSet scoreNegSet( motif, bgModel, negSet );
        scoreNegSet.calcLogOdds();
        std::vector<std::vector<float>> negAllScores = scoreNegSet.getMopsScores();
        std::vector<float> negScores;
        for( size_t i = 0; i < negSet.size(); i++ ){
            negScores.insert( std::end( negScores ),
                              std::begin( negAllScores[i] ),
                              std::end( negAllScores[i] ) );
        }

        // score positive sequence set
        // calculate p-values based on positive and negative scores
        ScoreSeqSet scorePosSet( motif, bgModel, posSet );
        scorePosSet.calcLogOdds();

        std::vector<std::vector<float>> posScores = scorePosSet.getMopsScores();
        scorePosSet.calcPvalues( posScores, negScores );

        scorePosSet.write( GScan::outputDirectory,
                           GScan::outputFileBasename + fileExtension,
                           GScan::pvalCutoff,
                           GScan::ss );

        delete motif;
    }

    if( bgModel ) delete bgModel;
    GScan::destruct();

    return 0;
}

