#include <time.h>		// time(), clock_t, clock, CLOCKS_PER_SEC
#include <stdio.h>

#include "Global.h"
#include "../shared/BackgroundModel.h"
#include "../shared/utils.h"
#include "MotifSet.h"
#include "EM.h"
#include "FDR.h"

int main( int nargs, char* args[] ){

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

    if( Global::debugMode ){
      Global::debug();
      Alphabet::debug();
    }

	fprintf( stderr, "\n" );
	fprintf( stderr, "************************\n" );
	fprintf( stderr, "*   Background Model   *\n" );
	fprintf( stderr, "************************\n" );
	BackgroundModel* bgModel = new BackgroundModel( *Global::negSequenceSet,
													Global::bgModelOrder,
													Global::bgModelAlpha,
													Global::interpolateBG );

	fprintf( stderr, "\n" );
	fprintf( stderr, "*********************\n" );
	fprintf( stderr, "*   Initial Motif   *\n" );
	fprintf( stderr, "*********************\n" );
	MotifSet motifs;

	// learn motifs
	fprintf( stderr, "\n" );
	fprintf( stderr, "*************************\n" );
	fprintf( stderr, "*   Learn Motif by EM   *\n" );
	fprintf( stderr, "*************************\n" );
	for( int N = 0; N < motifs.getN(); N++ ){
		Motif* motif = new Motif( *motifs.getMotifs()[N] );
		EM em( motif, bgModel );
		em.learnMotif();
		em.write();
	}

	// write motifs
	motifs.write();

	// evaluate motifs
	if( Global::FDR ){
		fprintf( stderr, "\n" );
		fprintf( stderr, "***********\n" );
		fprintf( stderr, "*   FDR   *\n" );
		fprintf( stderr, "***********\n" );
		for( int N = 0; N < motifs.getN(); N++ ){
			Motif* motif = new Motif( *motifs.getMotifs()[N] );
			FDR fdr( motif );
			fdr.evaluateMotif();
			fdr.write();
		}
	}

	fprintf( stderr, "\n" );
	fprintf( stderr, "******************\n" );
	fprintf( stderr, "*   Statistics   *\n" );
	fprintf( stderr, "******************\n" );
	std::cout << "Given alphabet type is " << Alphabet::getAlphabet();
	// for positive sequence set
	std::cout << "\nGiven positive sequence set is " << Global::posSequenceBasename
			<< "\n	"<< Global::posSequenceSet->getN() << " sequences. max.length: " <<
			Global::posSequenceSet->getMaxL() << ", min.length: " <<
			Global::posSequenceSet->getMinL() << "\n	base frequencies:";
	for( int i = 0; i < Alphabet::getSize(); i++ ){
		std::cout << ' ' << Global::posSequenceSet->getBaseFrequencies()[i]
		          << "(" << Alphabet::getAlphabet()[i] << ")";
	}
	// for negative sequence set
	std::cout << "\nGiven negative sequence set is " << Global::negSequenceBasename
			<< "\n	"<< Global::negSequenceSet->getN() << " sequences. max.length: "
			<< Global::negSequenceSet->getMaxL() << ", min.length: " <<
			Global::negSequenceSet->getMinL() << "\n	base frequencies:";
	for( int i = 0; i < Alphabet::getSize(); i++ )
		std::cout << ' ' << Global::negSequenceSet->getBaseFrequencies()[i]
				  << "(" << Alphabet::getAlphabet()[i] << ")";

	std::cout << "\nGiven initial model is " << Global::initialModelBasename;
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
