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
                                       GScan::posSequenceBasename );

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

    /**
     * Initialize the model
     */

    MotifSet motif_set( GScan::initialModelFilename,
                        GScan::addColumns.at(0),
                        GScan::addColumns.at(1),
                        GScan::initialModelTag,
                        GScan::posSequenceSet,
                        GScan::posSequenceSet->getBaseFrequencies(),
                        GScan::modelOrder,
                        GScan::modelAlpha );

	size_t motifNum = ( GScan::maxPWM > motif_set.getN() ) ? motif_set.getN() : GScan::maxPWM;

    /**
     * Sample negative sequence set based on s-mer frequencies
     */

    std::vector<Sequence*>  negset;
    size_t mFold = 10;
    // sample negative sequence set B1set based on s-mer frequencies
    // from positive training sequence set
    std::vector<std::unique_ptr<Sequence>> negSeqs;
    SeqGenerator generateNegSeqs( GScan::posSequenceSet->getSequences() );
    negSeqs = generateNegSeqs.arti_bgseqset( mFold );
    // convert unique_ptr to regular pointer
    for( size_t n = 0; n < negSeqs.size(); n++ ) {
        negset.push_back( negSeqs[n].release() );
    }

    for( size_t n = 0; n < motifNum; n++ ) {
        // deep copy each motif in the motif set
        Motif *motif = new Motif( *motif_set.getMotifs()[n] );

        // score negative sequence set
        ScoreSeqSet scoreNegSet( motif, bgModel, negset );
        scoreNegSet.calcLogOdds();
        std::vector<std::vector<float>> negAllScores = scoreNegSet.getMopsScores();
        std::vector<float> negScores;
        for( size_t n = 0; n < negset.size(); n++ ){
            negScores.insert( std::end( negScores ),
                              std::begin( negAllScores[n] ),
                              std::end( negAllScores[n] ) );
        }

        // score positive sequence set
        // calculate p-values based on positive and negative scores
        ScoreSeqSet scorePosSet( motif, bgModel, GScan::posSequenceSet->getSequences() );
        scorePosSet.calcLogOdds();
        std::vector<std::vector<float>> posScores = scorePosSet.getMopsScores();
        scorePosSet.calcPvalues( posScores, negScores );

        std::string fileExtension;
        if( GScan::initialModelTag == "PWM"){
            fileExtension = "_motif_" + std::to_string( n+1 );
        } else {
            fileExtension = GScan::fileExtension;
        }
        scorePosSet.write( GScan::outputDirectory,
                           GScan::posSequenceBasename + fileExtension,
                           GScan::pvalCutoff,
                           GScan::ss );

        delete motif;
    }
    if( bgModel ) delete bgModel;
    GScan::destruct();
    return 0;
}

