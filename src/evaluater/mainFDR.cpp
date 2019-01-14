//
// Created by wanwan on 03.09.17.
//

#include "../Global/Global.h"
#include "FDR.h"

int main( int nargs, char* args[] ){

    /**
     * initialization
     */

	Global G( nargs, args );

    /**
     * Build up the background model
     */
	BackgroundModel* bgModel;
    if( Global::bgModelFilename == NULL ) {
        bgModel = new BackgroundModel(Global::posSequenceSet->getSequences(),
                                      Global::outputFileBasename);
    } else {
        bgModel = new BackgroundModel( Global::bgModelFilename );
    }

    if( Global::saveInitialBaMMs ){
        // save background model
        bgModel->write(Global::outputDirectory, Global::outputFileBasename);
    }

    /**
     * Initialize the model
     */
    MotifSet motif_set( Global::initialModelFilename,
                        Global::addColumns.at(0),
                        Global::addColumns.at(1),
                        Global::initialModelTag,
                        Global::posSequenceSet );

    /**
     * Filter out short sequences
     */
    std::vector<Sequence*> posSet = Global::posSequenceSet->getSequences();
    size_t count = 0;
    std::vector<Sequence*>::iterator it = posSet.begin();
    while( it != posSet.end() ){
        if( (*it)->getL() < motif_set.getMaxW() ){
            if( Global::verbose ){
                std::cout << "Warning: remove the short sequence: "
                          << (*it)->getHeader() << std::endl;
            }
            posSet.erase(it);
            count++;
        } else {
            it++;
        }
    }
    if( count > 0 ){
        std::cout << "Note: " << count << " short sequences have been filtered out." << std::endl;
    }

    /**
     * Optional: subsample input sequence set when it is too large
     */
    size_t posN = posSet.size();
    if( Global::fixedPosN && Global::maxPosN < posN ){
        std::random_shuffle(posSet.begin(), posSet.end());
        for( size_t n = posN-Global::maxPosN; n > 0; n-- ){
            posSet.erase( posSet.begin() + Global::maxPosN + n - 1 );
        }
    }

    /**
     * Generate negative sequence set for cross-validation
     */
    std::vector<Sequence*>  negset;
    if( Global::B3 ){
        // take the given negative sequence set
        negset = Global::negSequenceSet->getSequences();

    } else {
        // generate negative sequence set based on s-mer frequencies
        // from positive training sequence set
        posN = posSet.size();   // update the size of positive sequences after filtering
        std::vector<std::unique_ptr<Sequence>> negSeqs;
        SeqGenerator negseq( posSet, NULL );
        if( !Global::fixedNegN and posN >= Global::negSeqNum ){
            negSeqs = negseq.sample_bgseqset_by_fold( Global::mFold );
            std::cout << Global::mFold << " x " << posN
                      << " background sequences are generated." << std::endl;
        } else if( !Global::fixedNegN and posN < Global::negSeqNum ){
            bool rest = Global::negSeqNum % posN;
            Global::mFold = Global::negSeqNum / posN + rest;
            negSeqs = negseq.sample_bgseqset_by_fold( Global::mFold );
            std::cout << Global::mFold << " x " << posN
                      << " background sequences are generated." << std::endl;
        } else {
            negSeqs = negseq.sample_bgseqset_by_num( Global::negSeqNum, Global::posSequenceSet->getMaxL() );
            std::cout << Global::negSeqNum << " (fixed) background sequences are generated." << std::endl;
        }

        // convert unique_ptr to regular pointer
        for (size_t n = 0; n < negSeqs.size(); n++) {
            negset.push_back(negSeqs[n].release());
            negSeqs[n].get_deleter();
        }
    }

    /**
     * Cross-validate the motif model
     */
    // Determine whether to optimize over motifs
    // or over optimization steps within a motif
    size_t mainLoopThreads;
    size_t perLoopThreads;
    if( Global::parallelOverMotif ){
        mainLoopThreads = Global::threads;
        perLoopThreads = 1;
    } else {
        mainLoopThreads = 1;
        perLoopThreads = Global::threads;
    }

#pragma omp parallel for num_threads( mainLoopThreads )

    for( size_t n = 0; n < motif_set.getN(); n++ ){

        Motif* motif = new Motif( *motif_set.getMotifs()[n] );

        FDR fdr( posSet, negset, motif, bgModel );

        fdr.evaluateMotif( perLoopThreads );

        if(Global::saveInitialBaMMs){
            // write out the foreground model
            motif->write( Global::outputDirectory,
                          Global::outputFileBasename + "_init_motif_" + std::to_string( n+1 ) );
        }

        std::string fileExtension;
        if( Global::initialModelTag == "PWM" ) {
            fileExtension = "_motif_" + std::to_string(n + 1);
        }

        fdr.write( Global::outputDirectory,
                   Global::outputFileBasename + fileExtension );
        if( motif )		delete motif;
    }

    // free memory
    if( !Global::B3 and negset.size() ) {
        for (size_t n = 0; n < negset.size(); n++) {
            delete negset[n];
        }
    }

    if( bgModel ) delete bgModel;

    return 0;
}