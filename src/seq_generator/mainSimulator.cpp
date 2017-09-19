//
// Created by wanwan on 18.09.17.
//

/*
 * This function is used for simulating sequences given different parameters.
 */

#include "GSimu.h"
#include "SeqGenerator.h"
#include "../init/MotifSet.h"
#include "../refinement/EM.h"

int main( int nargs, char* args[] ){

    /**
     * initialization
     */

    GSimu::init( nargs, args );

    if( GSimu::sampleBgset ){
        /**
         * Sample background sequence set based on s-mer frequencies from the input sequence set
         */
        SeqGenerator negseq( GSimu::sequenceSet->getSequences(), NULL, GSimu::sOrder );
        negseq.write( GSimu::outputDirectory,
                      GSimu::sequenceBasename + "_sampledNegset",
                      negseq.arti_negset( GSimu::mFold ) );
    }

    if( GSimu::maskSeqset or GSimu::embedSeqset ){
        /**
         * Build up the background model
         */
        BackgroundModel* bgModel = new BackgroundModel( GSimu::sequenceSet->getSequences(),
                                                        GSimu::bgModelOrder,
                                                        GSimu::bgModelAlpha );
        /**
         * Construct the motif set
         */
        MotifSet motif_set( GSimu::initialModelFilename,
                            GSimu::addColumns.at(0),
                            GSimu::addColumns.at(1),
                            GSimu::initialModelTag,
                            GSimu::sequenceSet,
                            GSimu::sequenceSet->getBaseFrequencies(),
                            GSimu::modelOrder,
                            GSimu::modelAlpha );

        for( size_t n = 0; n < motif_set.getN(); n++ ) {

            Motif* motif = new Motif( *motif_set.getMotifs()[n] );
            /**
             * optimize motif using EM
             */
            EM model( motif, bgModel, GSimu::sequenceSet->getSequences(), GSimu::q );
			model.optimize();

            if( GSimu::maskSeqset ){

                /**
                 * Mask the given motif from the input sequence set
                 */
                SeqGenerator seqset( GSimu::sequenceSet->getSequences(), motif, GSimu::modelOrder, GSimu::q );
                seqset.write( GSimu::outputDirectory,
                              GSimu::sequenceBasename + "_maskedPosset",
                              seqset.arti_negset_motif_masked( model.getR() ) );
            }

            if( GSimu::embedSeqset ){

                /**
                 * embed the given motif into the input sequence set
                 */
                SeqGenerator seqset( GSimu::sequenceSet->getSequences(), motif, GSimu::modelOrder, GSimu::q );
                seqset.write( GSimu::outputDirectory,
                              GSimu::sequenceBasename + "_embededSeqset",
                              seqset.arti_posset_motif_embedded( GSimu::mFold ) );
            }
            if( motif ) delete motif;
        }
        if( bgModel )   delete bgModel;
    }

    GSimu::destruct();

    return 0;
}