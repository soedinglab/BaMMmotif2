#include <iomanip>

#include "Global.h"
#include "BackgroundModel.h"
#include "MotifSet.h"
#include "EM.h"
#include "GibbsSampling.h"
#include "ScoreSeqSet.h"
#include "SeqGenerator.h"
#include "FDR.h"

int main( int nargs, char* args[] ){

	clock_t t0 = clock();

	std::cout << std::endl
              << "======================================" << std::endl
              << "=      Welcome to use BaMM!motif     =" << std::endl
              << "=                   Version 2.0      =" << std::endl
              << "=              by Soeding Group      =" << std::endl
              << "=  http://www.mpibpc.mpg.de/soeding  =" << std::endl
              << "======================================" << std::endl;

	// seed random number
	srand( 42 );
	Global::rngx.seed( 42 );

	// initialization
	Global::init( nargs, args );
    std::vector<Sequence*> posseqs = Global::posSequenceSet->getSequences();
    std::vector<Sequence*> negseqs = Global::negSequenceSet->getSequences();

	if( Global::verbose ){
		std::cout << std::endl
                  << "************************" << std::endl
                  << "*   Background Model   *" << std::endl
                  << "************************" << std::endl;
	}
	BackgroundModel* bgModel;
	if( !Global::bgModelGiven ){
		bgModel = new BackgroundModel( negseqs,
                                       Global::bgModelOrder,
                                       Global::bgModelAlpha,
                                       Global::interpolateBG,
                                       Global::posSequenceBasename );
	} else {
		bgModel = new BackgroundModel( Global::bgModelFilename );
	}

	// always save background model
	bgModel->write( Global::outputDirectory, Global::posSequenceBasename );

	if( Global::verbose ){
        std::cout << std::endl
                  << "***************************" << std::endl
                  << "*   Initial Motif Model   *" << std::endl
                  << "***************************" << std::endl;
	}

	MotifSet motif_set( Global::initialModelFilename,
                        Global::addColumns.at(0),
                        Global::addColumns.at(1),
                        Global::initialModelTag );

	size_t motifNum = ( Global::num > motif_set.getN() ) ? motif_set.getN() : Global::num;

	if( Global::verbose ){
        std::cout << std::endl
                  << "*********************" << std::endl
                  << "*   BaMM training   *" << std::endl
                  << "*********************" << std::endl;
	}

    for( size_t n = 0; n < motifNum; n++ ){
		// deep copy each motif in the motif set
		Motif* motif = new Motif( *motif_set.getMotifs()[n] );

		if( Global::saveInitialBaMMs ){
			// optional: save initial model
			motif->write( Global::outputDirectory,
                          Global::posSequenceBasename + "_init_motif_" + std::to_string( n+1 ) );
		}

		// optimize the model with either EM or Gibbs sampling
		if( Global::EM ){
			EM model( motif, bgModel, posseqs, Global::q );
			// learn motifs by EM
			model.optimize();
			// write model parameters on the disc
			if( Global::saveBaMMs ){
				model.write( Global::outputDirectory,
                             Global::posSequenceBasename + "_motif_" + std::to_string( n+1 ),
                             Global::ss );
			}
            // write out the learned model
            motif->write( Global::outputDirectory,
                          Global::posSequenceBasename + "_motif_" + std::to_string( n+1 ) );

		} else if ( Global::CGS ){
			GibbsSampling model( motif, bgModel, posseqs, Global::q );
			// learn motifs by collapsed Gibbs sampling
			model.optimize();
			// write model parameters on the disc
			if( Global::saveBaMMs ){
				model.write( Global::outputDirectory,
                             Global::posSequenceBasename + "_motif_" + std::to_string( n+1 ),
                             Global::ss );
			}
            // write out the learned model
            motif->write( Global::outputDirectory,
                          Global::posSequenceBasename + "_motif_" + std::to_string( n+1 ) );

		} else {

			std::cout << "Note: the model is not optimized!\n";
 //           motif->write( Global::outputDirectory, Global::posSequenceBasename + "_motif_ " + std::to_string( n+1 ) );
		}

		if( Global::scoreSeqset ){
			// score the model on sequence set
			ScoreSeqSet seq_set( motif, bgModel, posseqs );
			seq_set.score();
			seq_set.write( Global::outputDirectory,
						   Global::posSequenceBasename + "_motif_" + std::to_string( n+1 ),
						   Global::scoreCutoff,
                           Global::ss );
		}

		if( Global::generatePseudoSet ){

			// optimize motifs by EM
			EM model( motif, bgModel, posseqs, Global::q );
			model.optimize();

			// generate artificial sequence set with the learned motif embedded
			SeqGenerator artificial_seq_set( posseqs, motif, Global::sOrder );

			artificial_seq_set.write( Global::outputDirectory,
									  Global::posSequenceBasename + "_embed_motif_" + std::to_string( n+1 ),
									  artificial_seq_set.arti_posset_motif_embedded( Global::mFold ) );

		}

        if( Global::maskPosSequenceSet ){

            // learn motifs by EM
            EM model_opti( motif, bgModel, posseqs, Global::q );
            model_opti.optimize();

            // generate sequences from positive sequences with masked motif after optimization
            SeqGenerator artificial_set( posseqs, motif, Global::sOrder );
            artificial_set.write( Global::outputDirectory,
								  Global::posSequenceBasename + "_masked_motif_" + std::to_string( n+1 ),
                                  artificial_set.arti_negset_motif_masked( model_opti.getR() ) );
        }

        delete motif;
	}


	// evaluate motifs
	if( Global::FDR ){
		if( Global::verbose ){
            std::cout << std::endl
                      << "***********************" << std::endl
                      << "*   BaMM validation   *" << std::endl
                      << "***********************" << std::endl;
		}

        /**
         * Generate negative sequence set for cross-validation
         */

        std::vector<Sequence*>  negset;

        bool B1 = true;
        bool B2 = false;
        bool B3 = false;
        bool B3prime = false;

        if( B1 ) {
            // sample negative sequence set B1set based on s-mer frequencies
            // from positive training sequence set
            std::vector<std::unique_ptr<Sequence>> B1Seqs;
            SeqGenerator b1seqs( posseqs, NULL ,Global::sOrder );
            B1Seqs = b1seqs.arti_negset( Global::mFold );
            // convert unique_ptr to regular pointer
            for( size_t n = 0; n < B1Seqs.size(); n++ ) {
                negset.push_back( B1Seqs[n].release() );
            }
        }

        if( B2 ) {
            // sample negative sequence set B2set based on s-mer frequencies
            // from the given sampled negative sequence set
            std::vector<std::unique_ptr<Sequence>> B2Seqs;
            SeqGenerator b2seqs( negseqs, NULL ,Global::sOrder );
            B2Seqs = b2seqs.arti_negset( 1 );
            // convert unique_ptr to regular pointer
            for( size_t n = 0; n < B2Seqs.size(); n++ ) {
                negset.push_back( B2Seqs[n].release() );
            }
        }

        if( B3 ) {
            // take negative seqset as B3set
            negset = negseqs;
        }

        if( B3prime ) {
            // generate sequences from positive sequences with masked motif
            Motif *motif_opti = new Motif( *motif_set.getMotifs()[0] );
            SeqGenerator artificial_set( negseqs, motif_opti, Global::sOrder );
            EM model_opti( motif_opti, bgModel, posseqs, Global::q );
            // learn motifs by EM
            model_opti.optimize();
            std::vector<std::unique_ptr<Sequence>> B1SeqSetPrime;
            B1SeqSetPrime = artificial_set.arti_negset_motif_masked( model_opti.getR() );
            if ( motif_opti ) delete motif_opti;
            // draw B1setPrime from the positive sequence set with masked motifs
            for( size_t n = 0; n < B1SeqSetPrime.size(); n++ ) {
                negset.push_back(B1SeqSetPrime[n].release());
            }
        }

		// cross-validate the motif model
		for( size_t n = 0; n < motifNum; n++ ){
			Motif* motif = new Motif( *motif_set.getMotifs()[n] );
			FDR fdr( posseqs, negset, Global::q, motif, bgModel, Global::cvFold );
			fdr.evaluateMotif();
			fdr.write( Global::outputDirectory,
                       Global::posSequenceBasename + "_motif_" + std::to_string( n+1 ) );
			if( motif )		delete motif;
		}

	}

    std::cout << std::endl
              << "******************" << std::endl
              << "*   Statistics   *" << std::endl
              << "******************" << std::endl;
	Global::printStat();

	fprintf( stdout, "\n------ Runtime: %.2f seconds (%0.2f minutes) -------\n",
			( ( float )( clock() - t0 ) ) / CLOCKS_PER_SEC,
			( ( float )( clock() - t0 ) ) / ( CLOCKS_PER_SEC * 60.0f ) );

	// free memory
	if( bgModel ) delete bgModel;

	Global::destruct();

	return 0;
}
