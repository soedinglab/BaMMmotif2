#ifndef BACKGROUNDMODELSET_H_
#define BACKGROUNDMODELSET_H_

#include <fstream>	// e.g. std::ifstream
#include <vector>

#include <dirent.h>		// e.g. opendir
#include <string.h>		// e.g. strcmp
#include <sys/stat.h>

#include "BackgroundModel.h"
#include "SequenceSet.h"
#include "../refinement/utils.h"

class BackgroundModelSet{

public:

	// learn background init from sequences
	BackgroundModelSet( char* indir, char* extension, size_t order,
						std::vector<float> alpha, bool interpolate );
	// read in background init from files
	BackgroundModelSet( char* indir, char* extension );

	~BackgroundModelSet();

	std::vector<BackgroundModel*>& getBackgroundModels();
	size_t getN();

	// calculate log likelihoods for the sequence set
	// afterwards the background init contain log probabilities in v_
	std::vector<double> calculateLogLikelihoods( std::vector<Sequence*> seqs );

	// calculate posterior probabilities for the sequence set
	std::vector<double> calculatePosteriorProbabilities( std::vector<Sequence*>
															seqs );

	// calculate positional likelihoods for the sequence set
	// and write likelihoods to file
	void calculatePosLikelihoods( std::vector<Sequence*> seqs, char* odir );

	void print();
	void write( char* odir );

private:

	std::vector<BackgroundModel*>	backgroundModels_;
};

#endif /* BACKGROUNDMODELSET_H_ */
