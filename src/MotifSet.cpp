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

		length += Global::addColumns.at(0) + Global::addColumns.at(1);

		Motif* motif = new Motif( length );
		motif->initFromBindingSites( Global::bindingSitesFilename );
		motifs_.push_back( motif );
		N_ = 1;

	} else if( Global::PWMFilename != NULL ){

		// read file to calculate motif length
		// Motif* motif = new Motif( length )
		// motif.initFromPWM( file )
		// motifs.push_back( motif )

	} else if( Global::BaMMFilename != NULL ){

		// read file to calculate motif length
		// Motif* motif = new Motif( length )
		// motif.initFromBMM( file )
		// motifs.push_back( motif )
	}

//	if( Global::verbose ) print();
}

MotifSet::~MotifSet(){
//	motifs_.clear();
}

std::vector<Motif*> MotifSet::getMotifs(){
    return motifs_;
}

int MotifSet::getN(){
    return N_;
}

void MotifSet::print(){
	printf( " ____________________________\n"
			"|                            |\n"
			"| PROBABILITIES for MotifSet |\n"
			"|____________________________|\n\n" );

	int N = 0;
	for( int i = 0; i < N_; i++ ){
		printf( " ________________________________________\n"
				"|                                        |\n"
				"| INITIALIZED PROBABILITIES for Motif %d  |\n"
				"|________________________________________|\n\n", N+1 );
		motifs_[i]->print();
		N++;
	}
}
void MotifSet::write(){
	// before writing, remove the previous file
	std::string opath = std::string( Global::outputDirectory )  + '/'
			+ std::string( Global::posSequenceBasename ) + ".condsInit";
	remove( opath.c_str() );					// avoid appending if the file already exists
	for( int i = 0; i < N_; i++ )
		motifs_[i]->write();
}

