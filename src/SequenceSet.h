#ifndef SEQUENCESET_H_
#define SEQUENCESET_H_

#include <fstream>	// e.g. std::ifstream
#include <limits>	// e.g. std::numeric_limits
#include <sstream>	// e.g. std::ostringstream
#include <vector>

#include <math.h>	// e.g. logf

#include "Alphabet.h"
#include "Sequence.h"
#include "utils.h"

class SequenceSet{

public:

	SequenceSet( std::string sequenceFilepath,
			bool singleStrand = false,
			std::string intensityFilepath = "",
			float q = 0.9f );
	~SequenceSet();

	std::string				getSequenceFilepath();
	std::string				getIntensityFilepath();
	std::vector<Sequence*> 	getSequences();
	size_t 					getN();
	size_t 					getMinL();
	size_t					getMaxL();
	float					getQ();
	float* 					getBaseFrequencies();

	void					print();			// print sequences

private:

	std::string				sequenceFilepath_;	// path to FASTA file
	std::string				intensityFilepath_;	// path to intensity file

	std::vector<Sequence*>	sequences_;			// sequences
	size_t 					N_;					// number of sequences
	size_t 					minL_;				// length of the shortest sequence
	size_t 					maxL_;				// length of the longest sequence
	float					q_;					// the fraction of sequences with
												// motif
	float*	 				baseFrequencies_;	// kmer frequencies

	std::vector<size_t>		Y_;					// contains 1 at position 0
												// and the number of oligomers y for increasing order k at positions k+1
												// e.g.
												// alphabet size_ = 4: Y_ = 4^0 4^1 4^2 ... 4^15 < std::numeric_limits<int>::max()
												// limits the length of oligomers to 15 (and the order to 14)

	int 					readFASTA( bool singleStrand = false );	// read in FASTA file
	int 					readIntensities();	// read in intensity file
};

#endif /* SEQUENCESET_H_ */
