/*
 * SequenceSet.h
 *
 *  Created on: Apr 1, 2016
 *      Author: wanwan
 */

#ifndef SEQUENCESET_H_
#define SEQUENCESET_H_

#include <stdlib.h>

#include "Sequence.h"

class SequenceSet {

public:

	SequenceSet( char* sequenceFilename, char* intensityFilename = NULL ){
		/*
		 * read in sequences
		 * calculate mono-nucleotide frequencies
		 * read in intensities (optional)
		 */
	}
	~SequenceSet();

	char* 		getFilename();				// get input sequence filename
	Sequence* getSequences();				// get sequences
	int 			getN();								// get number of sequences
	int 			getMinL();						// get min. length of sequences
	int 			getMaxL();						// get max. length of sequences
	float* 		getBaseFrequencies();	// get mono-nucleotide frequencies

private:

	char* 		filename_;						// input sequence filename
	Sequence* sequences_;						// sequences
	int 			N_;										// number of sequences
	int 			minL_;								// min. length of sequences
	int 			maxL_;								// max. length of sequences
	float* 		baseFrequencies_;			// mono-nucleotide frequencies

	int 			readFASTA_();					// read in FASTA file
	int 			readIntensities_();		// read in intensity file
};

#endif /* SEQUENCESET_H_ */
