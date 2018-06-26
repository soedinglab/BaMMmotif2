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
                      GSimu::sequenceBasename + "_bgset",
                      negseq.sample_bgseqset_by_fold(GSimu::mFold) );
    } else {
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
                            bgModel->getV(),
                            GSimu::bgModelOrder,
                            GSimu::modelOrder,
                            GSimu::modelAlpha );

        for( size_t n = 0; n < motif_set.getN(); n++ ) {

            Motif *motif = new Motif(*motif_set.getMotifs()[n]);

            if (GSimu::maskSeqset) {
                /**
                 * optimize motif using EM
                 */
                EM model(motif, bgModel, GSimu::sequenceSet->getSequences(), GSimu::q);
                model.EStep();
                /**
                 * Mask the given motif from the input sequence set
                 */
                SeqGenerator seq_generator(GSimu::sequenceSet->getSequences(),
                                           motif,
                                           GSimu::modelOrder,
                                           GSimu::q);
                seq_generator.write(GSimu::outputDirectory,
                                    GSimu::sequenceBasename + "_motif_" + std::to_string(n + 1) + "_masked",
                                    seq_generator.seqset_with_motif_masked(model.getR()));
            } else if (GSimu::embedSeqset) {
                /**
                 * embed the given motifs into the input sequence set
                 */
                SeqGenerator seq_generator(GSimu::sequenceSet->getSequences(),
                                           motif,
                                           GSimu::modelOrder,
                                           GSimu::q);
                seq_generator.write(GSimu::outputDirectory,
                                    GSimu::sequenceBasename + "_motif_" + std::to_string(n + 1) + "_embedded",
                                    seq_generator.arti_posset_motif_embedded(GSimu::at));
            } else {
                std::cout << "No artificial sequence set is generated. Please check your input options."
                          << std::endl;
            }
            if (motif) delete motif;

        }
        if( bgModel )   delete bgModel;
    }

    GSimu::destruct();

    return 0;
}