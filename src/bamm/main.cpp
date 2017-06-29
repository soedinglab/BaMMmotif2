#include <iomanip>
#include <time.h>		// time(), clock_t, clock, CLOCKS_PER_SEC
#include <stdio.h>

#include "Global.h"
#include "../shared/BackgroundModel.h"
#include "../shared/utils.h"
#include "MotifSet.h"
#include "ModelLearning.h"
#include "ScoreSeqSet.h"
#include "SeqGenerator.h"
#include "FDR.h"

int main( int nargs, char* args[] ){

	clock_t t0 = clock();

	// seed random number
	srand( 42 );
	Global::rngx.seed( 42 );

	fprintf( stderr, "\n" );
	fprintf( stderr, "======================================\n" );
	fprintf( stderr, "=      Welcome to use BaMMmotif      =\n" );
	fprintf( stderr, "=                   Version 1.0      =\n" );
	fprintf( stderr, "=              by Soeding Group      =\n" );
	fprintf( stderr, "=  http://www.mpibpc.mpg.de/soeding  =\n" );
	fprintf( stderr, "======================================\n" );

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
		bgModel = new BackgroundModel( *Global::bgSequenceSet,
										Global::bgModelOrder,
										Global::bgModelAlpha,
										Global::interpolateBG );
	} else {
		bgModel = new BackgroundModel( Global::bgModelFilename );
	}

	if( Global::saveBgModel ){
		bgModel->write( Global::outputDirectory );
	}

	if( Global::verbose ){
		fprintf( stderr, "\n" );
		fprintf( stderr, "**************************\n" );
		fprintf( stderr, "*   Initial Motif Model  *\n" );
		fprintf( stderr, "**************************\n" );
	}

	MotifSet motifs;

	size_t motifNum = ( Global::num > motifs.getN() ) ? motifs.getN() : Global::num;

	if( Global::verbose ){
		fprintf( stderr, "\n" );
		fprintf( stderr, "********************\n" );
		fprintf( stderr, "*   BaMM training  *\n" );
		fprintf( stderr, "********************\n" );
	}

	for( size_t n = 0; n < motifNum; n++ ){
		// deep copy each motif in the motif set
		Motif* motif = new Motif( *motifs.getMotifs()[n] );

		// train the model with either EM or Gibbs sampling
		ModelLearning model( motif, bgModel );

		if( Global::EM ){

			// learn motifs by EM
			model.EM();

		} else if ( Global::CGS ){

			// learn motifs by collapsed Gibbs sampling
			model.GibbsSampling();

		} else {

			std::cout << "Model is not optimized!\n";

			//exit( -1 );	// allow to do FDR on Initial Model i.e. for PWMs or precalculated BaMMs

		}

		if( Global::generatePseudoSet ){

			// optimize motifs by EM
			model.EM();

			// generate pesudo positive sequence set based on the learned motif
			SeqGenerator seqset( Global::posSequenceSet->getSequences(), model.getMotif() );

			seqset.sample_pseudo_seqset( Global::mFold );

			seqset.write_pseudoset();

		}

		// write model parameters on the disc
		if( Global::saveBaMMs ){
			model.write( n+1 );
		}

		// write out the learned model
		motif->write( n+1 );

		if( Global::scoreSeqset ){
			// score the model on sequence set
			ScoreSeqSet seqset( motif, bgModel, Global::posSequenceSet->getSequences() );
			seqset.score();
			seqset.write( Global::outputDirectory, Global::scoreCutoff );
		}

		delete motif;
	}

	// evaluate motifs
	if( Global::FDR ){
		if( Global::verbose ){
			fprintf( stderr, "\n" );
			fprintf( stderr, "*********************\n" );
			fprintf( stderr, "*   Validate BaMM   *\n" );
			fprintf( stderr, "*********************\n" );
		}
		for( size_t n = 0; n < motifNum; n++ ){
			Motif* motif = new Motif( *motifs.getMotifs()[n] );
			FDR fdr( motif, Global::cvFold );
			fdr.evaluateMotif( n );
			fdr.write( n );
			delete motif;
		}
	}

	fprintf( stderr, "\n" );
	fprintf( stderr, "******************\n" );
	fprintf( stderr, "*   Statistics   *\n" );
	fprintf( stderr, "******************\n" );
	Global::printStat();

	fprintf( stdout, "\n-------------- Runtime: %.2f seconds (%0.2f minutes) --------------\n",
			( ( float )( clock() - t0 ) ) / CLOCKS_PER_SEC, ( ( float )( clock() - t0 ) ) / ( CLOCKS_PER_SEC * 60.0f ) );

	// free memory
	if( bgModel ) delete bgModel;
	Global::destruct();

	return 0;
}
