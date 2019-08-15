//
// Created by wanwan on 13.08.19.
//

#include "../Global/Global.h"
#include "../init/MotifSet.h"
#include "KmerPredicter.h"

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
        if( (*it)->getL() < motif_set.getMaxW() || (*it)->getL() < Global::kmerLength ){
            std::cout << "Warning: remove the short sequence: " << (*it)->getHeader()
                      << std::endl;
            posSet.erase(it);
        } else {
            it++;
        }
    }

    /**
     * Count k-mers
     */
    KmerPredicter Kd( Global::kmerLength );
    Kd.countKmer( posSet );

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

        /**
         * Score k-mers
         */
        Kd.scoreKmer( motif, bgModel );
        Kd.writeKmerCounts( Global::outputDirectory,
                            Global::outputFileBasename + fileExtension );

        delete motif;
    }

    delete bgModel;

    return 0;
}

