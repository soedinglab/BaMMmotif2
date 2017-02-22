#ifndef BACKGROUNDMODELSET_H_
#define BACKGROUNDMODELSET_H_

#include <fstream>	// e.g. std::ifstream
#include <vector>

#include <dirent.h>		// e.g. opendir
#include <string.h>		// e.g. strcmp
#include <sys/stat.h>

#include "BackgroundModel.h"
#include "SequenceSet.h"
#include "utils.h"

class BackgroundModelSet{

public:

	// learn background models from sequences
	BackgroundModelSet( char* inputDirectory, char* extension, int order, std::vector<float> alpha, bool interpolate );
	// read in background models from files
	BackgroundModelSet( char* inputDirectory, char* extension );

	~BackgroundModelSet();

	std::vector<BackgroundModel*>& getBackgroundModels();
	size_t getN();

	// calculate log likelihoods for the sequence set
	// afterwards the background models contain log probabilities in v_
	std::vector<double> calculateLogLikelihoods( SequenceSet& sequenceSet );

	// calculate posterior probabilities for the sequence set
	std::vector<double> calculatePosteriorProbabilities( SequenceSet& sequenceSet );

	// calculate positional likelihoods for the sequence set
	// and write likelihoods to file
	void calculatePosLikelihoods( SequenceSet& sequenceSet, char* outputDirectory );

	void print();
	void write( char* dir );

private:

	std::vector<BackgroundModel*>	backgroundModels_;
};

#endif /* BACKGROUNDMODELSET_H_ */
