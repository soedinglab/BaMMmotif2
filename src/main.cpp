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
#include "EM.h"
#include "FDR.h"

int main( int nargs, char* args[] ){

	long timestamp = time( NULL );

	// initialization
	Global::init( nargs, args );

	MotifSet motifs;	// here need to check if the motif length exceeds the min. posSeq length!

	BackgroundModel bg;
	std::vector< std::vector<int> > folds;
	bg.init( folds );	// not sure about the input folds parameter

/*
	// learn motifs
	for( std::list<Motif*>::const_iterator iter = motifs.getMotifs.begin(); iter != motifs.getMotifs.end(); iter++ ){
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
*/
	fprintf( stderr, "It seems that Global destructor doesn't work. \n" );
//	Global::destruct();

	fprintf( stdout, "\nRuntime: %ld seconds (%0.2f minutes)\n", time( NULL )-timestamp,
			( float )( time( NULL )-timestamp )/60.0f );

	return 0;
}
