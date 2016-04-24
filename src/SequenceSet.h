/*
 * SequenceSet.h
 *
 *  Created on: Apr 1, 2016
 *      Author: wanwan
 */

#ifndef SEQUENCESET_H_
#define SEQUENCESET_H_

#include <vector>			//std::vector

#include <stdlib.h>			// malloc, calloc, free, rand, exit

#include "Sequence.h"

class SequenceSet {

public:

	SequenceSet( char* sequenceFilepath, char* intensityFilepath = NULL );

	~SequenceSet();

	char* 					getSequenceFilepath();	            // get input sequence filename
	char* 					getIntensityFilepath();             // get input intensity filename
	std::vector<Sequence*> 	getSequences();			            // get sequences
	unsigned int 			getN();					            // get number of sequences
	unsigned int 			getMinL();				            // get min. length of sequences
	unsigned int			getMaxL();				            // get max. length of sequences
	float* 					getBaseFrequencies();	            // get mono-nucleotide frequencies

private:

	char*					sequenceFilepath_;		            // input sequence filename
	char*					intensityFilepath_;		            // input intensity filename
	std::vector<Sequence*>	sequences_;							// sequences
	unsigned int 			N_;						            // number of sequences
	unsigned int 			minL_;					            // min. length of sequences
	unsigned int 			maxL_;					            // max. length of sequences
	float* 					baseFrequencies_;		            // mono-nucleotide frequencies

	int 					readFASTA();			            // read in FASTA file
	int 					readIntensities();		            // read in intensity file
};

#endif /* SEQUENCESET_H_ */
