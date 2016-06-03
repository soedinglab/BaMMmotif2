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
	fprintf( stderr, "======================================\n" );
	fprintf( stderr, "=      Welcome to use BaMMmotif      =\n" );
	fprintf( stderr, "=                   Version 1.0      =\n" );
	fprintf( stderr, "=              by Soeding Group      =\n" );
	fprintf( stderr, "=  http://www.mpibpc.mpg.de/soeding  =\n" );
	fprintf( stderr, "======================================\n" );

	// set up the precision for all float or double values
	std::cout << std::setprecision(4) << std::endl;

	// initialization
	Global::init( nargs, args );

	fprintf( stderr, "\n" );
	fprintf( stderr, "************************\n" );
	fprintf( stderr, "*   Background Model   *\n" );
	fprintf( stderr, "************************\n" );
	BackgroundModel bgModel;
	bgModel.init( std::vector<int> () );

	fprintf( stderr, "\n" );
	fprintf( stderr, "*********************\n" );
	fprintf( stderr, "*   Initial Motif   *\n" );
	fprintf( stderr, "*********************\n" );

	MotifSet motifs;

	// learn motifs
	fprintf( stderr, "\n" );
	fprintf( stderr, "**********\n" );
	fprintf( stderr, "*   EM   *\n" );
	fprintf( stderr, "**********\n" );
	for( int i = 0; i < motifs.getN(); i++ ){
		std::cout << "N_ = " << motifs.getN() << std::endl;
		EM em( motifs.getMotifs()[i], bgModel, std::vector<int> () );
		em.learnMotif();
		std::cout <<" Marker 1 " << std::endl;
		em.write();
		std::cout <<" Marker 2 " << std::endl;

	}
	exit( -1 );
	std::cout <<" Marker 3 " << std::endl;


	// write motifs
	motifs.write();
	/*
	// evaluate motifs
	fprintf( stderr, "\n" );
	fprintf( stderr, "***********\n" );
	fprintf( stderr, "*   FDR   *\n" );
	fprintf( stderr, "***********\n" );
	for( std::list<Motif*>::const_iterator iter = motifs.getMotifs().begin(); iter != motifs.getMotifs().end(); iter++ ){
		FDR fdr( *iter );
		fdr.evaluateMotif();
		fdr.write();
	}
*/
//	Global::destruct();

	fprintf( stderr, "\n" );
	fprintf( stderr, "******************\n" );
	fprintf( stderr, "*   Statistics   *\n" );
	fprintf( stderr, "******************\n" );
	std::cout << "Chosen positive sequence set is " << Global::posSequenceFilename << std::endl;
	if(Global::negSequenceFilename )
		std::cout << "Chosen negative sequence set is " << Global::negSequenceFilename << std::endl;
	std::cout << "Chosen alphabet type is " << Alphabet::getAlphabet() << std::endl;
	std::cout << "Folds for FDR estimation: " << Global::cvFold << std::endl;

	fprintf( stdout, "\nRuntime: %ld seconds (%0.2f minutes)\n", time( NULL )-timestamp,
			( float )( time( NULL )-timestamp )/60.0f );

	return 0;
}
