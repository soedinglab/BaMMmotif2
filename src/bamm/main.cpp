#include <time.h>		// time()
#include <stdio.h>

#include "Global.h"
#include "../shared/BackgroundModel.h"
#include "../shared/utils.h"
#include "MotifSet.h"
#include "EM.h"
#include "FDR.h"

int main( int nargs, char* args[] ){

	long timestamp = time( NULL );

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
													Global::bgModelAlpha );
	bgModel->write( Global::outputDirectory );

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
		EM em( motifs.getMotifs()[N], bgModel );
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
			FDR fdr( motifs.getMotifs()[N] );
			fdr.evaluateMotif();
			fdr.write();
		}
	}

	fprintf( stderr, "\n" );
	fprintf( stderr, "******************\n" );
	fprintf( stderr, "*   Statistics   *\n" );
	fprintf( stderr, "******************\n" );
	std::cout << "Given alphabet type is " << Alphabet::getAlphabet();
	std::cout << "\nGiven positive sequence set is " << Global::posSequenceBasename
			<< "\n	"<< Global::posSequenceSet->getN() << " sequences. max.length: " <<
			Global::posSequenceSet->getMaxL() << ", min.length: " <<
			Global::posSequenceSet->getMinL() << "\n	base frequencies:";
	for( int i = 0; i < Alphabet::getSize(); i++ ){
		std::cout << ' ' << Global::posSequenceSet->getBaseFrequencies()[i]
		          << "(" << Alphabet::getAlphabet()[i] << ")";
	}
	if( Global::negSequenceFilename ){
		std::cout << "\nGiven negative sequence set is " << Global::negSequenceBasename
				<< "\n	"<< Global::negSequenceSet->getN() << " sequences. max.length: "
				<< Global::negSequenceSet->getMaxL() << ", min.length: " <<
				Global::negSequenceSet->getMinL() << "\n	base frequencies:";
		for( int i = 0; i < Alphabet::getSize(); i++ )
			std::cout << ' ' << Global::negSequenceSet->getBaseFrequencies()[i]
			          << "(" << Alphabet::getAlphabet()[i] << ")";
	} else {
		std::cout << "\nNone negative sequence set is given";
	}
	std::cout << "\nGiven initial model is " << Global::initialModelBasename;
	if( Global::FDR ){
		std::cout << "\nGiven folds for FDR estimation: " << Global::cvFold;
	}
	if( Global::setSlow ){
		std::cout << "\n***** This is a slow EM version. *****";
	} else {
		std::cout << "\n***** This is a fast EM version. *****";
	}

	fprintf( stdout, "\n-------------- Runtime: %ld seconds (%0.2f minutes) --------------\n",
			time( NULL )-timestamp, ( float )( time( NULL )-timestamp )/60.0f );

	// free memory
	if( bgModel ) delete bgModel;
	Global::destruct();

	return 0;
}
