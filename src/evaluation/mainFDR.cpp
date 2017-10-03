//
// Created by wanwan on 03.09.17.
//


#include "GFdr.h"
#include "FDR.h"

int main( int nargs, char* args[] ){

    /**
     * initialization
     */

	GFdr::init( nargs, args );

    /**
     * Build up the background model
     */

	BackgroundModel* bgModel = new BackgroundModel( GFdr::negSequenceSet->getSequences(),
                                                    GFdr::bgModelOrder,
                                                    GFdr::bgModelAlpha,
                                                    GFdr::interpolateBG,
                                                    GFdr::posSequenceBasename );
    /**
     * Initialize the model
     */

    MotifSet motif_set( GFdr::initialModelFilename,
                        GFdr::addColumns.at(0),
                        GFdr::addColumns.at(1),
                        GFdr::initialModelTag,
                        GFdr::posSequenceSet,
                        GFdr::negSequenceSet->getBaseFrequencies(),
                        GFdr::modelOrder,
                        GFdr::modelAlpha );

	size_t motifNum = ( GFdr::num > motif_set.getN() ) ? motif_set.getN() : GFdr::num;

    /**
     * Generate negative sequence set for cross-validation
     */

    std::vector<Sequence*>  negset;
    // generate negative sequence set based on s-mer frequencies
    // from positive training sequence set
    std::vector<std::unique_ptr<Sequence>> negSeqs;
    SeqGenerator negseq( GFdr::posSequenceSet->getSequences(), NULL, GFdr::sOrder );
    negSeqs = negseq.arti_bgseqset(GFdr::mFold);
    // convert unique_ptr to regular pointer
    for( size_t n = 0; n < negSeqs.size(); n++ ) {
        negset.push_back( negSeqs[n].release() );
    }

    /**
     * Cross-validate the motif model
     */

    for( size_t n = 0; n < motifNum; n++ ){
        Motif* motif = new Motif( *motif_set.getMotifs()[n] );
        FDR fdr( GFdr::posSequenceSet->getSequences(), negset, GFdr::q,
                 motif, bgModel,
                 GFdr::cvFold, GFdr::mops, GFdr::zoops,
                 true, GFdr::savePvalues, GFdr::saveLogOdds );

        fdr.evaluateMotif( GFdr::EM, GFdr::CGS );
        fdr.write( GFdr::outputDirectory,
                   GFdr::posSequenceBasename + "_motif_" + std::to_string( n+1 ) );
        if( motif )		delete motif;
    }

    GFdr::destruct();

    return 0;
}