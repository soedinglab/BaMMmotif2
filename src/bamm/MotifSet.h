/*
 * MotifSet.h
 *
 *  Created on: Apr 1, 2016
 *      Author: wanwan
 */

#ifndef MOTIFSET_H_
#define MOTIFSET_H_

#include "Motif.h"

class MotifSet {

public:

	MotifSet();
	~MotifSet();

	std::vector<Motif*> getMotifs();		// get motifs
	int			 		getN();				// get number of motifs

	void 				print();			// print motifs to console
	void 				write();			// write motifs to files (basename.bmm)

private:

	std::vector<Motif*> motifs_;			// motifs
	int          		N_ = 0;				// number of motifs
};

#endif /* MOTIFSET_H_ */
