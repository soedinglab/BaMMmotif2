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
	std::cout << "Destructor for MotifSet class works fine. \n";
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

	for( int i = 0; i < N_; i++ ){
		printf( " ________________________________________\n"
				"|                                        |\n"
				"| INITIALIZED PROBABILITIES for Motif %d  |\n"
				"|________________________________________|\n\n", i+1 );
		motifs_[i]->print();
	}
}
void MotifSet::write(){

	/**
	 * save all the motifs learned by BaMM in one flat flie:
	 * posSequenceBasename.motifs: list conditional probabilities for each motif after EM training
	 */

	std::string opath = std::string( Global::outputDirectory )  + '/'
			+ std::string( Global::posSequenceBasename ) + ".motifs";
	std::ofstream ofile( opath.c_str() );

	int i, j, k, y, W;
	float*** v_motif;

	for( i = 0; i < N_; i++ ){
		ofile << "Motif " << i+1 <<":" << std::endl;
		W =  motifs_[i]->getW();
		v_motif =  motifs_[i]->getV();
		for( j = 0; j < W; j++ ){
			for( k = 0; k < Global::modelOrder+1; k++ ){
				for( y = 0; y < Global::powA[k+1]; y++ )
					ofile << std::scientific << std::setprecision(8) << v_motif[k][y][j] << '\t';
				ofile << std::endl;
			}
			ofile << std::endl;
		}
	}
}
