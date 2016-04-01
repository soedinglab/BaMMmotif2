/*
 * main.cpp
 *
 *  Created on: Apr 1, 2016
 *      Author: wanwan
 */

#include "Global.h"
#include "MotifSet.h"
#include "EM.h"
#include "FDR.h"

int main( int nargs, char *args[] ){

	long timestamp = time( NULL );

	// learn motifs

	// write motifs

	// evaluate motifs

	fprintf( stdout, "\nRuntime: %ld seconds (%0.2f minutes)\n", time( NULL )-timestamp,
			( float )( time( NULL )-timestamp )/60.0f );
	return 0;
}



