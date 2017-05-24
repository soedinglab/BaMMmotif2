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

	// seed random number
	srand( 42 );

	clock_t t0 = clock();

	fprintf( stderr, "\n" );
	fprintf( stderr, "======================================\n" );
	fprintf( stderr, "=      Welcome to use BaMMmotif      =\n" );
	fprintf( stderr, "=                   Version 1.0      =\n" );
	fprintf( stderr, "=              by Soeding Group      =\n" );
	fprintf( stderr, "=  http://www.mpibpc.mpg.de/soeding  =\n" );
	fprintf( stderr, "======================================\n" );

	// initialization
	Global::init( nargs, args );

	Global::rngx.seed( 42 );

	fprintf( stderr, "\n" );
	fprintf( stderr, "************************\n" );
	fprintf( stderr, "*   Background Model   *\n" );
	fprintf( stderr, "************************\n" );
	BackgroundModel* bgModel;

	if( !Global::bgModelGiven ){
	bgModel = new BackgroundModel( *Global::negSequenceSet,
													Global::bgModelOrder,
													Global::bgModelAlpha,
													Global::interpolateBG );
	} else {
		bgModel = new BackgroundModel( Global::bgModelFilename );
	}

	if( Global::saveBgModel ){
		bgModel->write( Global::outputDirectory );
	}

	fprintf( stderr, "\n" );
	fprintf( stderr, "*********************\n" );
	fprintf( stderr, "*   Initial Motif   *\n" );
	fprintf( stderr, "*********************\n" );

	MotifSet motifs;

	int motifNum = ( Global::num > motifs.getN() ) ? Global::num : motifs.getN();

	for( int n = 0; n < motifNum; n++ ){
		// initialize the model
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

			seqset.sample_pseudo_seqset();

			seqset.write( seqset.sample_pseudo_seqset() );

		}

		// write model parameters on the disc
		if( Global::saveBaMMs ){
			model.write( n );
		}

		// write out the learned model
		motif->write( n );

		if ( Global::scoreSeqset){
			// score the model on sequence set
			ScoreSeqSet seqset( motif, bgModel, Global::posSequenceSet->getSequences() );
			seqset.score();
			seqset.write( n, Global::scoreCutoff );
		}

		delete motif;
	}

	// evaluate motifs
	if( Global::FDR ){
		fprintf( stderr, "\n" );
		fprintf( stderr, "***********\n" );
		fprintf( stderr, "*   FDR   *\n" );
		fprintf( stderr, "***********\n" );
		for( int n = 0; n < motifNum; n++ ){
			Motif* motif = new Motif( *motifs.getMotifs()[n] );
			FDR fdr( motif );
			fdr.evaluateMotif( n );
			fdr.write( n );
			delete motif;
		}
	}

	fprintf( stderr, "\n" );
	fprintf( stderr, "******************\n" );
	fprintf( stderr, "*   Statistics   *\n" );
	fprintf( stderr, "******************\n" );
	std::cout << "Given alphabet type is " << Alphabet::getAlphabet();
	std::cout << "\nGiven initial model is " << Global::initialModelBasename;
	// for positive sequence set
	std::cout << "\nGiven positive sequence set is " << Global::posSequenceBasename
			<< "\n	"<< Global::posSequenceSet->getN() << " sequences. max.length: " <<
			Global::posSequenceSet->getMaxL() << ", min.length: " <<
			Global::posSequenceSet->getMinL() << "\n	base frequencies:";
	for( int i = 0; i < Alphabet::getSize(); i++ ){
		std::cout << ' ' << Global::posSequenceSet->getBaseFrequencies()[i]
		          << "(" << Alphabet::getAlphabet()[i] << ")";
	}
	std::cout << "\n	" << Global::q * 100 << "% of the sequences contain the optimized motif.";
	// for negative sequence set
	if( Global::negSeqGiven ){
		std::cout << "\nGiven negative sequence set is " << Global::negSequenceBasename
				<< "\n	"<< Global::negSequenceSet->getN() << " sequences. max.length: "
				<< Global::negSequenceSet->getMaxL() << ", min.length: " <<
				Global::negSequenceSet->getMinL() << "\n	base frequencies:";
		for( int i = 0; i < Alphabet::getSize(); i++ )
			std::cout << ' ' << Global::negSequenceSet->getBaseFrequencies()[i]
					  << "(" << Alphabet::getAlphabet()[i] << ")";
	}

	if( Global::FDR ){
		std::cout << "\nGiven folds for FDR estimation: " << Global::cvFold;
	}

	fprintf( stdout, "\n-------------- Runtime: %.2f seconds (%0.2f minutes) --------------\n",
			( ( float )( clock() - t0 ) ) / CLOCKS_PER_SEC, ( ( float )( clock() - t0 ) ) / ( CLOCKS_PER_SEC * 60.0f ) );

	// free memory
	if( bgModel ) delete bgModel;
	Global::destruct();

	return 0;
}
