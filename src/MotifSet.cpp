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

	} else if( Global::bindingSitesFilename ){

		std::cout << "Reading in the binding sites file: " << Global::bindingSitesFilename << std::endl;

		std::ifstream file;
		file.open( Global::bindingSitesFilename, std::ifstream::in );
		std::string seq;
		int length;

		if( !file.good() ){
			std::cout << "Error: Cannot open bindingSitesFile sequence file: "
					<< Global::bindingSitesFilename << std::endl;
			exit( -1 );
		} else {
			getline( file, seq );					// get length of the first sequence
			length = seq.length();
		}

		Motif* motif = new Motif( length );

		std::cout << "Initializing motif from bindingSitesFile..." << std::endl;

		motif->initFromBindingSites( Global::bindingSitesFilename );

		std::cout << "Motif from bindingSitesFile has been initialized." << std::endl;

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

	print();
	write();
}

MotifSet::~MotifSet(){

}

std::list<Motif*> MotifSet::getMotifs(){
    return motifs_;
}

int MotifSet::getN(){
    return N_;
}

void MotifSet::print(){
	for( int n = 0; n < N_; n++ ){

	}
}
void MotifSet::write(){

}

