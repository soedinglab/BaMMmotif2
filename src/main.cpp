/*
 * main.cpp
 *
 *  Created on: Apr 1, 2016
 *      Author: wanwan
 */

#include "MotifSet.h"
#include "EM.h"
#include "FDR.h"

#include <time.h>
#include <stdio.h>
#include <list>

int main( int nargs, char *args[] ){

	long timestamp = time( NULL );

	// initialization
	Global::init( nargs, args );

	MotifSet motifs;

	BackgroundModel bg;
	bg.init();

	// learn motifs
	for( std::list<Motif*>::const_iterator iter=motifs.getMotifs.begin(); iter != motifs.getMotifs.end(); iter++ ){
		EM em = new EM( *iter, bg );
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

	Global::destruct();

	fprintf( stdout, "\nRuntime: %ld seconds (%0.2f minutes)\n", time( NULL )-timestamp,
			( float )( time( NULL )-timestamp )/60.0f );
	return 0;
}



