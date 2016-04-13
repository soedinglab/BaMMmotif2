/*
 * SequenceSet.h
 *
 *  Created on: Apr 1, 2016
 *      Author: wanwan
 */

#ifndef SEQUENCESET_H_
#define SEQUENCESET_H_

#include <stdlib.h>			/* malloc, calloc, free, rand, exit */

#include "Sequence.h"

class SequenceSet {

public:

	SequenceSet( char* sequenceFilepath, char* intensityFilepath = NULL ){
	/*
	 * read in sequences
	 * calculate mono-nucleotide frequencies
	 * read in intensities (optional)
	*/
	    // read in FASTA file
	    readFASTA( sequenceFilepath );

	    // if intensity file is given, read in intensity file
	    if ( intensityFilepath != NULL )
	        readIntensities( intensityFilepath );
	}
	~SequenceSet();

	char* 			getFilename(){ return filename_; };				    // get input sequence filename
	Sequence* 		getSequences(){ return sequences_; };				// get sequences
	unsigned int 	getN(){ return N_; };								// get number of sequences
	unsigned int 	getMinL(){ return minL_; };							// get min. length of sequences
	unsigned int	getMaxL(){ return maxL_; };							// get max. length of sequences
	float* 			getBaseFrequencies(){ return baseFrequencies_; };	// get mono-nucleotide frequencies

private:

	char* 			filename_;					// input sequence filename
	Sequence* 		sequences_;					// sequences
	unsigned int 	N_;							// number of sequences
	unsigned int 	minL_;						// min. length of sequences
	unsigned int 	maxL_;						// max. length of sequences
	float* 			baseFrequencies_;			// mono-nucleotide frequencies

	int 			readFASTA( char* sequenceFilepath );				// read in FASTA file
	int 			readIntensities( char* intensityFilepath );		    // read in intensity file
};

#endif /* SEQUENCESET_H_ */
