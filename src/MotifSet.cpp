/*
 * MotifSet.cpp
 *
 *  Created on: Apr 18, 2016
 *      Author: administrator
 */

#include <fstream>		// std::fstream

#include "MotifSet.h"

MotifSet::MotifSet(){

	if( Global::BaMMpatternFilename != NULL ){

	    // scan file and conduct for IUPAC pattern p
		// * Motif* motif = new Motif( length(p) )
		// * motif.initFromIUPACPattern( p )
		// motifs.push_back( motif )

	} else if( Global::bindingSitesFilename != NULL ){

		std::ifstream file( Global::bindingSitesFilename );
		std::string seq;
		int length;

		if( !file ){
			fprintf( stderr, "Cannot open positive sequence file: %s\n",
					Global::bindingSitesFilename );
			exit( -1 );
		} else {
			getline( file, seq );					// get length of the first sequence
			length = seq.length();
		}

		Motif* motif = new Motif( length );

		motif->initFromBindingSites( Global::bindingSitesFilename );

		motifs_.push_back( motif );
		N_ = 1;

		delete motif;

	} else if( Global::PWMFilename != NULL ){

		// read file to calculate motif length
		// Motif* motif = new Motif( length )
		// motif.initFromPWM( file )
		// motifs.push_back( motif )

	} else if( Global::BMMFilename != NULL ){

		// read file to calculate motif length
		// Motif* motif = new Motif( length )
		// motif.initFromBMM( file )
		// motifs.push_back( motif )
	}

}

MotifSet::~MotifSet(){

}

std::list<Motif*> MotifSet::getMotifs(){
    return motifs_;
}

int MotifSet::getN(){
    return N_;
}

void MotifSet::print(){}
void MotifSet::write(){}

