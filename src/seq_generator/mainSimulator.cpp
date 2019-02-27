//
// Created by wanwan on 18.09.17.
//

/*
 * This function is used for simulating sequences given different parameters.
 */

#include "../Global/Global.h"
#include "SeqGenerator.h"
#include "../init/MotifSet.h"
#include "../refiner/EM.h"

int main( int nargs, char* args[] ){

    /**
     * initialization
     */
    Global G( nargs, args );

    if( Global::sampleBgset ){
        /**
         * Sample background sequence set based on s-mer frequencies
         * from the input sequence set
         */
        SeqGenerator negseq( Global::posSequenceSet->getSequences(), NULL);
        negseq.write( Global::outputDirectory,
                      Global::posSequenceBasename + "_bgset",
                      negseq.sample_bgseqset_by_fold(Global::mFold) );
    } else {
        /**
         * Build up the background model
         */
        BackgroundModel* bgModel = new BackgroundModel( Global::negSequenceSet->getSequences(),
                                                        Global::outputFileBasename );
        /**
         * Construct the motif set
         */
        MotifSet motif_set( Global::initialModelFilename,
                            Global::addColumns.at(0),
                            Global::addColumns.at(1),
                            Global::initialModelTag,
                            Global::posSequenceSet );

        for( size_t n = 0; n < motif_set.getN(); n++ ) {

            Motif *motif = new Motif( *motif_set.getMotifs()[n] );

            if ( Global::maskSeqset ) {
                /**
                 * optimize motif using EM
                 */
                EM model( motif, bgModel, Global::posSequenceSet->getSequences() );
                model.EStep();
                /**
                 * Mask the given motif from the input sequence set
                 */
                SeqGenerator seq_generator( Global::posSequenceSet->getSequences(), motif );
                seq_generator.write(Global::outputDirectory,
                                    Global::posSequenceBasename + "_motif_" + std::to_string(n+1) + "_masked",
                                    seq_generator.seqset_with_motif_masked( model.getR() ));
            } else if ( Global::embedSeqset ) {
                /**
                 * embed the given motifs into the input sequence set
                 */
                SeqGenerator seq_generator(Global::posSequenceSet->getSequences(), motif);
                seq_generator.write(Global::outputDirectory,
                                    Global::posSequenceBasename + "_motif_" + std::to_string(n+1) + "_embedded",
                                    seq_generator.arti_posset_motif_embedded( Global::at ));
            } else {
                std::cout << "No artificial sequence set is generated. "
                          << "Please check your input options." << std::endl;
            }
            if (motif) delete motif;

        }
        if( bgModel )   delete bgModel;
    }

    return 0;
}