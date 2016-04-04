/*
 * SequenceSet.h
 *
 *  Created on: Apr 1, 2016
 *      Author: wanwan
 */

#ifndef SEQUENCESET_H_
#define SEQUENCESET_H_

#include "Sequence.h"

#include <stdlib.h>

class SequenceSet {

public:

	SequenceSet( char* sequenceFilename, char* intensityFilename=NULL ){
		// read in sequences
		// calculate mono-nucleotide frequencies
		// read in intensities (optional)
	}
	~SequenceSet();

	char* 		getFilename();				// get input sequence filename
	Sequence* 	getSequences();				// get sequences
	int 		getN();						// get number of sequences
	int 		getMinL();					// get min. length of sequences
	int 		getMaxL();					// get max. length of sequences
	float* 		getBaseFrequencies();		// get mono-nucleotide frequencies

private:

	char* 		filename;					// input sequence filename
	Sequence* 	sequences;					// sequences
	int 		N;							// number of sequences
	int 		minL;						// min. length of sequences
	int 		maxL;						// max. length of sequences
	float* 		baseFrequencies;			// mono-nucleotide frequencies

	int 		readFASTA();				// read in FASTA file
	int 		readIntensities();			// read in intensity file
};

#endif /* SEQUENCESET_H_ */
