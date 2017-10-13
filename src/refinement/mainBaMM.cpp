#include <iomanip>

#include "Global.h"
#include "EM.h"
#include "GibbsSampling.h"
#include "../evaluation/FDR.h"

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

	if( Global::verbose ){
		std::cout << std::endl
                  << "************************" << std::endl
                  << "*   Background Model   *" << std::endl
                  << "************************" << std::endl;
	}

	BackgroundModel* bgModel;
	if( !Global::bgModelGiven ){
		bgModel = new BackgroundModel( Global::posSequenceSet->getSequences(),
                                       Global::bgModelOrder,
                                       Global::bgModelAlpha,
                                       Global::interpolateBG,
                                       Global::posSequenceBasename );
        // always save background model
        bgModel->write( Global::outputDirectory, Global::posSequenceBasename );

    } else {
		bgModel = new BackgroundModel( Global::bgModelFilename );
		if( Global::savebgModel){
			bgModel->write( Global::outputDirectory, Global::posSequenceBasename );
		}
	}

	if( Global::verbose ){
        std::cout << std::endl
                  << "***************************" << std::endl
                  << "*   Initial Motif Model   *" << std::endl
                  << "***************************" << std::endl;
	}

	MotifSet motif_set( Global::initialModelFilename,
                        Global::addColumns.at(0),
                        Global::addColumns.at(1),
                        Global::initialModelTag,
                        Global::posSequenceSet,
                        Global::posSequenceSet->getBaseFrequencies(),
                        Global::modelOrder,
                        Global::modelAlpha );

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
			EM model( motif, bgModel, Global::posSequenceSet->getSequences(),
                      Global::q, Global::optimizeQ, Global::verbose );
			// learn motifs by EM
			if( !Global::advanceEM ) {
                model.optimize();
            } else {
                model.advance();
            }

            // write model parameters on the disc
			if( Global::saveBaMMs ){
				model.write( Global::outputDirectory,
                             Global::posSequenceBasename + "_motif_" + std::to_string( n+1 ),
                             Global::ss );
			}

            // print out the optimized q for checking:
            std::cout << "optimized q = " << model.getQ() << std::endl;

		} else if ( Global::CGS ){
			GibbsSampling model( motif, bgModel, Global::posSequenceSet->getSequences(), Global::q, !Global::noQSampling );
			// learn motifs by collapsed Gibbs sampling
			model.optimize();
			// write model parameters on the disc
			if( Global::saveBaMMs ){
				model.write( Global::outputDirectory,
                             Global::posSequenceBasename + "_motif_" + std::to_string( n+1 ),
                             Global::ss );
			}

            // print out the optimized q for checking:
            std::cout << "optimized q = " << model.getQ() << std::endl;

		} else {

			std::cout << "Note: the model is not optimized!\n";
		}
        // write out the learned model
        motif->write( Global::outputDirectory,
                      Global::posSequenceBasename + "_motif_" + std::to_string( n+1 ) );

        if( Global::scoreSeqset ){
            // score the model on sequence set
            ScoreSeqSet seq_set( motif, bgModel, Global::posSequenceSet->getSequences() );

            seq_set.score();
            seq_set.write( Global::outputDirectory,
                           Global::posSequenceBasename + "_motif_" + std::to_string( n+1 ),
                           Global::scoreCutoff,
                           Global::ss );
        }

        if ( Global::bammSearch ){
        			fprintf( stderr, "\n" );
            		fprintf( stderr, "**********************\n" );
        			fprintf( stderr, "*   Score Sequences  *\n" );
        			fprintf( stderr, "**********************\n" );

        			BackgroundModel* bg = bgModel;
        			// 0. Depending on motif input, and learning, define bgModel
        			// use bgModel generated from input sequences when prediction is turned on

        			if( ! Global::EM and ! Global::CGS ){

        				// use provided bgModelFile if initialized with bamm format
        				if( Global::BaMMFilename != NULL ){
        					if( Global::bgModelFile == NULL ){
        						std::cout << "No background Model File provided for initial search motif!\n";
        						exit(-1);
        					}
        					bg = new BackgroundModel( Global::bgModelFile );
        				}
        				else{
        					std::cout << "not optimized with BammFile\n";
        					// use bgModel generated when reading in PWM File
        					if( Global::PWMFilename != NULL ){
        						//std::cout << "Init with PWM instead \n";
        						bg = new BackgroundModel( Global::PWMFilename , 0, 1);
        						// this means that also the global motif order needs to be adjusted;
        						Global::modelOrder = 0;
        					}
        					else{
        						std::cout << "No background Model found for provided initial search motif!\n";
        						exit(-1);
        					}
        				}

        			}

            		std::cout << "Seqgenerator, neg seqs... ";
        			SeqGenerator neg_seqs( Global::posSequenceSet->getSequences(), model.getMotif() );
        			std::cout << ".. sample seqset... ";

        			//neg_seqs.sample_negative_seqset(  );
        			neg_seqs.arti_bgseqset(FOLD);

        			std::cout << "DONE!\n";

        			// generate and score negative sequence set
        			std::cout << "generate negative sequence set... ";
        			ScoreSeqSet neg_seqset( model.getMotif(), bg, neg_seqs.getSeqs());
        			std::cout << "... and score it ...";
        			neg_seqset.score();
        			std::cout << "DONE!\n";

        			// <------this is already done in the ScoreSeqset part so probably I can scip all that

        			// score positive sequence set
        			std::cout << "generate positive sequence set... ";
        			ScoreSeqSet pos_seqset( model.getMotif(), bg, Global::posSequenceSet->getSequences());
        			std::cout << "... and score it ...";
        			pos_seqset.score();
        			std::cout << "DONE!\n";

        			// collapse and rank negative scores
        			std::vector<float> neg_scores = neg_seqset.getScoreAll();
        			std::sort( neg_scores.begin(), neg_scores.end(), std::greater<float>() );

        			fprintf( stderr, "\n" );
            		fprintf( stderr, "***************************\n" );
            		fprintf( stderr, "*   Evaluate Occurrences  *\n" );
        			fprintf( stderr, "***************************\n" );

        			// calculate p- and e-values for positive scores based on negative scores
        			pos_seqset.calcPvalues(neg_scores);

        			fprintf( stderr, "\n" );
            		fprintf( stderr, "***************************\n" );
            		fprintf( stderr, "*   print PValues         *\n" );
        			fprintf( stderr, "***************************\n" );

        			// print out scores and p-/e-values larger than p-values based cutoff
        			pos_seqset.writePvalues( n, Global:: pvalCutoff );
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

        if( Global::B2 ) {
            // sample negative sequence set B2set based on s-mer frequencies
            // from the given sampled negative sequence set
            std::vector<std::unique_ptr<Sequence>> B2Seqs;
            SeqGenerator b2seqs( Global::negSequenceSet->getSequences() );
            B2Seqs = b2seqs.arti_bgseqset(1);
            // convert unique_ptr to regular pointer
            for( size_t n = 0; n < B2Seqs.size(); n++ ) {
                negset.push_back( B2Seqs[n].release() );
            }
        } else if( Global::B3 ) {
            // take negative seqset as B3set
            negset = Global::negSequenceSet->getSequences();
        } else if( Global::B3prime ) {
            // generate sequences from positive sequences with masked motif
            Motif *motif_opti = new Motif( *motif_set.getMotifs()[0] );
            SeqGenerator artificial_set( Global::negSequenceSet->getSequences(), motif_opti );
            EM model_opti( motif_opti, bgModel, Global::posSequenceSet->getSequences(), Global::q );
            // learn motifs by EM
            model_opti.optimize();
            std::vector<std::unique_ptr<Sequence>> B1SeqSetPrime;
            B1SeqSetPrime = artificial_set.seqset_with_motif_masked(model_opti.getR());
            if ( motif_opti ) delete motif_opti;
            // draw B1setPrime from the positive sequence set with masked motifs
            for( size_t n = 0; n < B1SeqSetPrime.size(); n++ ) {
                negset.push_back(B1SeqSetPrime[n].release());
            }
        } else {
            // sample negative sequence set B1set based on s-mer frequencies
            // from positive training sequence set
            std::vector<std::unique_ptr<Sequence>> B1Seqs;
            SeqGenerator b1seqs( Global::posSequenceSet->getSequences() );
            B1Seqs = b1seqs.arti_bgseqset(Global::mFold);
            // convert unique_ptr to regular pointer
            for( size_t n = 0; n < B1Seqs.size(); n++ ) {
                negset.push_back( B1Seqs[n].release() );
            }
        }

        /**
         * cross-validate the motif model
         */
		for( size_t n = 0; n < motifNum; n++ ){
			Motif* motif = new Motif( *motif_set.getMotifs()[n] );
			FDR fdr( Global::posSequenceSet->getSequences(), negset,
                     Global::q, motif, bgModel, Global::cvFold,
                     Global::mops, Global::zoops,
                     Global::savePRs, Global::savePvalues, Global::saveLogOdds );
			fdr.evaluateMotif(Global::EM, Global::CGS, Global::optimizeQ, Global::advanceEM);
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
