#ifndef VIRUSHOSTINTERACTIONS_H_
#define VIRUSHOSTINTERACTIONS_H_

#include <algorithm>	// e.g. std::sort
#include <limits>		// e.g. std::numeric_limits
#include <vector>

#include <dirent.h>		// e.g. opendir

#include "../shared/BackgroundModelSet.h"
#include "../shared/SequenceSet.h"
#include "../shared/utils.h"

class VirusHostInteractions{

public:

	VirusHostInteractions( char* inputDirectoryBaMMs, char* extensionBaMMs );
	~VirusHostInteractions();

	void predict( char* inputDirectorySeqs, char* extensionSeqs );

	void print();
	void write( char* outputDirectory );

private:

	void				aggregatePosteriors();										// aggregate M x N matrix of posterior probabilities
	std::vector<double> calculatePosteriors( std::vector<double> llikelihoods );	// calculate posterior probabilities from log likelihoods
	double				score( BackgroundModel& bamm, SequenceSet& sequenceSet );	// calculate log likelihood

	BackgroundModelSet* bamms_;

	std::vector<std::string> bammNames_;				// names of M background models
	std::vector<std::string> sequenceSetNames_;			// basenames of N sequence set files

	std::vector< std::vector<double> > posteriors_;		// M x N matrix of posterior probabilities

	// aggregate statistics
	std::vector<size_t> bestBammIndices_;				// indices of background models with highest posterior per sequence set
	std::vector<double> bestBammPosteriors_;			// highest posterior per sequence set
	std::vector<size_t> secondBestBammIndices_;			// indices of background models with second highest posterior per sequence set
	std::vector<double> secondBestBammPosteriors_;		// second highest posterior per sequence set
};

#endif /* VIRUSHOSTINTERACTIONS_H_ */
