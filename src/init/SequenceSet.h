#ifndef SEQUENCESET_H_
#define SEQUENCESET_H_

#include <fstream>	// e.g. std::ifstream
#include <limits>	// e.g. std::numeric_limits
#include <sstream>	// e.g. std::ostringstream
#include <vector>

#include <math.h>	// e.g. logf

#include "Alphabet.h"
#include "Sequence.h"
#include "../Global/utils.h"

class SequenceSet{

public:

	SequenceSet( std::string sequenceFilepath,
			bool singleStrand = false,
			std::string intensityFilepath = "" );

	~SequenceSet();

	std::string				getSequenceFilepath();
	std::string				getIntensityFilepath();
	std::vector<Sequence*> 	getSequences();
    std::vector<std::string> sequence2string();
	size_t 					getMinL();
	size_t					getMaxL();
    size_t                  getBaseSum();
	float* 					getBaseFrequencies();
    bool                    isSingleStranded();
	void					print();			// print sequences

private:

	std::string				sequenceFilepath_;	// path to FASTA file
	std::string				intensityFilepath_;	// path to intensity file

	std::vector<Sequence*>	sequences_;			// sequences
	size_t 					minL_;				// length of the shortest sequence
	size_t 					maxL_;				// length of the longest sequence
    size_t                  baseSum_;           // sum of all the bases
	float*	 				baseFrequencies_;	// kmer frequencies
    bool                    isSingleStranded_;  // flag for searching on single strand

	int 					readFASTA();        // read in FASTA file
	int 					readIntensities();	// read in intensity file
};

#endif /* SEQUENCESET_H_ */
