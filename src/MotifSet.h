/*
 * MotifSet.h
 *
 *  Created on: Apr 1, 2016
 *      Author: wanwan
 */

#ifndef MOTIFSET_H_
#define MOTIFSET_H_

#include "Global.h"
#include "Motif.h"

class MotifSet {

public:

	MotifSet(){

		if( Global::GIMMEpatternFilename != NULL ){

		    // scan file and conduct for IUPAC pattern p
			// * Motif* motif = new Motif( length(p) )
			// * motif.initFromIUPACPattern( p )
			// motifs.push_back( motif )

		} else if( Global::bindingSitesFilename != NULL ){

			// read file to calculate motif length
			// Motif* motif = new Motif( length )
			// motif.initFromBindingSites( file )
			// motifs.push_back( motif )

		} else if( Global::PWMFilename != NULL ){

			// read file to calculate motif length
			// Motif* motif = new Motif( length )
			// motif.initFromPWM( file )
			// motifs.push_back( motif )

		} else if( Global::iIMMFilename != NULL ){

			// read file to calculate motif length
			// Motif* motif = new Motif( length )
			// motif.initFromPWM( file )
			// motifs.push_back( motif )
		}
	}
	~MotifSet();

	list<Motif*> getMotifs();	// get motifs
	int			 getN();		// get number of motifs

	void print(); 				// print motifs to console
	void write();				// write motifs to files (basename.iimm)

private:

	list<Motif*> motifs;		// motifs
	int          N;				// number of motifs
};

#endif /* MOTIFSET_H_ */
