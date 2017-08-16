#include <iomanip>

#include "Global.h"
#include "BackgroundModel.h"
#include "MotifSet.h"
#include "ModelLearning.h"
#include "ScoreSeqSet.h"
#include "SeqGenerator.h"
#include "FDR.h"

int main( int nargs, char* args[] ){

	clock_t t0 = clock();

	fprintf( stderr, "\n" );
	fprintf( stderr, "======================================\n" );
	fprintf( stderr, "=      Welcome to use BaMM!motif     =\n" );
	fprintf( stderr, "=                   Version 2.0      =\n" );
	fprintf( stderr, "=              by Soeding Group      =\n" );
	fprintf( stderr, "=  http://www.mpibpc.mpg.de/soeding  =\n" );
	fprintf( stderr, "======================================\n" );

	// seed random number
	srand( 42 );
	Global::rngx.seed( 42 );

	// initialization
	Global::init( nargs, args );

	if( Global::verbose ){
		fprintf( stderr, "\n" );
		fprintf( stderr, "************************\n" );
		fprintf( stderr, "*   Background Model   *\n" );
		fprintf( stderr, "************************\n" );
	}
	BackgroundModel* bgModel;
	if( !Global::bgModelGiven ){
		bgModel = new BackgroundModel( Global::negSequenceSet->getSequences(),
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
		fprintf( stderr, "\n" );
		fprintf( stderr, "**************************\n" );
		fprintf( stderr, "*   Initial Motif Model  *\n" );
		fprintf( stderr, "**************************\n" );
	}

	MotifSet motif_set( Global::initialModelFilename,
					Global::addColumns.at(0),
					Global::addColumns.at(1),
					Global::initialModelTag );

	size_t motifNum = ( Global::num > motif_set.getN() ) ?
						motif_set.getN() : Global::num;

	if( Global::verbose ){
		fprintf( stderr, "\n" );
		fprintf( stderr, "********************\n" );
		fprintf( stderr, "*   BaMM training  *\n" );
		fprintf( stderr, "********************\n" );
	}

	for( size_t n = 0; n < motifNum; n++ ){
		// deep copy each motif in the motif set
		Motif* motif = new Motif( *motif_set.getMotifs()[n] );

		// train the model with either EM or Gibbs sampling
		ModelLearning model( motif,
							bgModel,
							Global::posSequenceSet->getSequences(),
							Global::q );

		if( Global::saveInitialBaMMs ){
			// optional: save initial model
			motif->write( Global::outputDirectory,
							Global::posSequenceBasename, 0 );
		}

		if( Global::EM ){
			// learn motifs by EM
			model.EM();

		} else if ( Global::CGS ){
			// learn motifs by collapsed Gibbs sampling
			model.GibbsSampling();

		} else {

			std::cout << "Note: the model is not optimized!\n";

		}

		// write model parameters on the disc
		if( Global::saveBaMMs ){
			model.write( Global::outputDirectory,
							Global::posSequenceBasename, n+1, Global::ss );
		}

		// write out the learned model
		// motif->write( n+1 );
		model.getMotif()->write( Global::outputDirectory,
									Global::posSequenceBasename, n+1 );

		if( Global::scoreSeqset ){
			// score the model on sequence set
			ScoreSeqSet seq_set( motif, bgModel,
								Global::posSequenceSet->getSequences() );
			seq_set.score();
			seq_set.write( Global::outputDirectory, Global::posSequenceBasename,
					n+1, Global::scoreCutoff, Global::ss );
		}

		if( Global::generatePseudoSet ){

			// optimize motifs by EM
			model.EM();

			// generate artificial sequence set with the learned motif embedded
			SeqGenerator seq_set( Global::posSequenceSet->getSequences(),
									motif,
									Global::sOrder );

			seq_set.write( Global::outputDirectory,
							Global::posSequenceBasename,
							n+1,
							seq_set.arti_posset_motif_embedded( Global::mFold ) );

		}

		delete motif;
	}

	// evaluate motifs
	if( Global::FDR ){
		if( Global::verbose ){
			fprintf( stderr, "\n" );
			fprintf( stderr, "*********************\n" );
			fprintf( stderr, "*  BaMM validation  *\n" );
			fprintf( stderr, "*********************\n" );
		}

		// cross-validate the motif model
		for( size_t n = 0; n < motifNum; n++ ){
			Motif* motif = new Motif( *motif_set.getMotifs()[n] );
			FDR fdr( Global::posSequenceSet->getSequences(),
					Global::negSequenceSet->getSequences(),
					Global::q,
					motif,
					Global::cvFold );
			fdr.evaluateMotif();
			fdr.write( Global::outputDirectory,
					Global::posSequenceBasename, n+1 );
			if( motif )		delete motif;
		}

	}

	fprintf( stderr, "\n" );
	fprintf( stderr, "******************\n" );
	fprintf( stderr, "*   Statistics   *\n" );
	fprintf( stderr, "******************\n" );
	Global::printStat();

	fprintf( stdout, "\n------ Runtime: %.2f seconds (%0.2f minutes) -------\n",
			( ( float )( clock() - t0 ) ) / CLOCKS_PER_SEC,
			( ( float )( clock() - t0 ) ) / ( CLOCKS_PER_SEC * 60.0f ) );

	// free memory
	if( bgModel ) delete bgModel;
	Global::destruct();

	return 0;
}
