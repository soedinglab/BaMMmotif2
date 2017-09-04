#ifndef BACKGROUNDMODELSETSCORE_H_
#define BACKGROUNDMODELSETSCORE_H_

#include <vector>

#include <dirent.h>		// e.g. opendir

#include "BackgroundModelSet.h"
#include "SequenceSet.h"
#include "../refinement/utils.h"

class BackgroundModelSetScore{

public:

	BackgroundModelSetScore( char* inputDirectoryBaMMs, char* extensionBaMMs );
	~BackgroundModelSetScore();

	// calculate posterior probabilities
	void predict( char* inputDirectorySeqs, char* extensionSeqs );
	// calculate positional likelihoods and write likelihoods to file
	void score( char* inputDirectorySeqs, char* extensionSeqs, char* outputDirectory );

	void print();
	void write( char* outputDirectory );

private:

	void aggregate();								// aggregate M x N matrix of
													// posterior probabilities

	BackgroundModelSet* 		bamms_;
	std::vector<std::string> 	bammNames_;			// names of M bg init
	std::vector<std::string> 	sequenceSetNames_;	// basenames of N sequence
													// set files

	std::vector<std::vector<double>> posteriors_;	// M x N matrix of posterior
													// probabilities

	// aggregate statistics
	std::vector<size_t> bestBammIndices_;			// indices of background init
													// with highest posterior per sequence set
	std::vector<double> bestBammPosteriors_;		// highest posterior per sequence set
	std::vector<size_t> secondBestBammIndices_;		// indices of background init with
													// second highest posterior per sequence set
	std::vector<double> secondBestBammPosteriors_;	// second highest posterior per sequence set
};

#endif /* BACKGROUNDMODELSETSCORE_H_ */
