/*
 * main.cpp
 *
 *  Created on: Apr 1, 2016
 *      Author: wanwan
 */

#include <list>			// list, list::begin, list::end
#include <vector>

#include <time.h>		// time()
#include <stdio.h>

#include "Global.h"
#include "MotifSet.h"
#include "BackgroundModel.h"
#include "EM.h"
#include "FDR.h"

int main( int nargs, char* args[] ){

	long timestamp = time( NULL );

	fprintf( stderr, "\n" );
	fprintf( stderr, "=================================\n" );
	fprintf( stderr, "=   Welcome to use BaMM MOTIF   =\n" );
	fprintf( stderr, "=================================\n" );

	// initialization
	Global::init( nargs, args );

	fprintf( stderr, "\n" );
	fprintf( stderr, "******************\n" );
	fprintf( stderr, "*   Statistics   *\n" );
	fprintf( stderr, "******************\n" );
	std::cout << "Chosen positive sequence set is " << Global::posSequenceFilename << std::endl;
	//std::cout << "Chosen negative sequence set is " << Global::negSequenceFilename << std::endl;
	std::cout << "Chosen alphabet type is " << Alphabet::getAlphabet() << std::endl;
	//std::cout << "Chosen binding site file is:" << Global::bindingSitesFilename << std::endl;
	std::cout << "Folds for FDR estimation: " << Global::nFolds << std::endl;

	BackgroundModel bgModel;
	bgModel.init();

	fprintf( stderr, "\n" );
	fprintf( stderr, "*********************\n" );
	fprintf( stderr, "*   Initial Motif   *\n" );
	fprintf( stderr, "*********************\n" );

	MotifSet motifs;
	/*
	// learn motifs
	for( std::list<Motif*>::iterator iter = motifs.getMotifs().begin(); iter != motifs.getMotifs().end(); iter++ ){
		EM em = new EM( *iter, bgModel );
		em.learnMotif();
		em.write();
	}

	// write motifs
	motifs.write();

	// evaluate motifs
	for( std::list<Motif*>::const_iterator iter=motifs.getMotifs.begin(); iter != motifs.getMotifs.end(); iter++ ){
		FDR fdr = new FDR( *iter );
		fdr.evaluateMotif();
		fdr.write();
	}
*/

//	Global::destruct();

	fprintf( stdout, "\nRuntime: %ld seconds (%0.2f minutes)\n", time( NULL )-timestamp,
			( float )( time( NULL )-timestamp )/60.0f );

	return 0;
}
