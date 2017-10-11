//
// Created by wanwan on 04.09.17.
// This function is used for scanning motif occurrences in a sequence set
//

#include "GScan.h"
#include "ScoreSeqSet.h"
#include "../init/MotifSet.h"
#include "../refinement/Global.h"

int main( int nargs, char* args[] ) {

    /**
     * initialization
     */

	GScan::init( nargs, args );

    /**
     * Build up the background model
     */

	BackgroundModel* bgModel = new BackgroundModel( GScan::negSequenceSet->getSequences(),
                                                    GScan::bgModelOrder,
                                                    GScan::bgModelAlpha,
                                                    GScan::interpolateBG,
                                                    GScan::posSequenceBasename );
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

	size_t motifNum = ( GScan::num > motif_set.getN() ) ? motif_set.getN() : GScan::num;

    for( size_t n = 0; n < motifNum; n++ ) {
        // deep copy each motif in the motif set
        Motif *motif = new Motif( *motif_set.getMotifs()[n] );
        // score the model on sequence set
        ScoreSeqSet seq_set( motif, bgModel, GScan::posSequenceSet->getSequences() );

        seq_set.score();
        seq_set.write(GScan::outputDirectory,
                      GScan::posSequenceBasename + "_motif_" + std::to_string( n + 1 ),
                      GScan::cutoff,
                      GScan::ss);
        delete motif;
    }

    GScan::destruct();
    return 0;
}

