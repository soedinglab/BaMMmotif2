#ifndef BACKGROUNDMODEL_H_
#define BACKGROUNDMODEL_H_

#include <fstream>	// e.g. std::ofstream
#include <iostream>	// e.g. std::cerr
#include <limits>	// e.g. std::numeric_limits
#include <numeric>	// e.g. std::iota
#include <string>
#include <vector>
#include <sys/stat.h>

#include <math.h>	// e.g. logf

#include "Alphabet.h"
#include "SequenceSet.h"
#include "utils.h"

class BackgroundModel{

public:

	BackgroundModel( SequenceSet& sequenceSet,
			         int order,
			         std::vector<float> alpha,
			         bool interpolate = true,
			         std::vector< std::vector<int> > foldIndices = std::vector< std::vector<int> >(),
			         std::vector<int> folds = std::vector<int>() );

	BackgroundModel( std::string filePath );

	BackgroundModel( char* , int K, float A );

	~BackgroundModel();

	std::string getName();
	int			getOrder();
	float**		getV();

	void 		expV();
	void 		logV();

	bool		vIsLog();

	// calculate the log likelihood for the sequence set
	// afterwards v_ contains log probabilities
	double calculateLogLikelihood( SequenceSet& sequenceSet,
			                       std::vector<std::vector<int>> foldIndices = std::vector<std::vector<int>>(),
			                       std::vector<int> folds = std::vector<int>() );

	// calculate positional likelihoods for the sequence set
	// and write likelihoods to file
	void calculatePosLikelihoods( SequenceSet& sequenceSet,
			                      char* outputDirectory,
			                      std::vector<std::vector<int>> foldIndices = std::vector<std::vector<int>>(),
			                      std::vector<int> folds = std::vector<int>() );


	void 	print();
    // afterwards v_ contains exp probabilities
	void 	write( char* dir );

private:

	// calculate probabilities from conditional probabilities
	void	calculateProbabilities( float** p );
	// calculate conditional probabilities from counts
	void 	calculateV();

	std::string			name_;					// basename of sequence set file

	int**				n_;						// oligomer counts
	float** 			v_;						// oligomer conditional probabilities

	bool				vIsLog_ = false;		// v_ contains log probabilities

	int					K_;						// order
	std::vector<float>	A_;						// order-specific alphas

	bool				interpolate_ = true;	// calculate prior probabilities from lower-order probabilities
												// instead of background frequencies of mononucleotides

	std::vector<int>	Y_;						// contains 1 at position 0
												// and the number of oligomers y for increasing order k (from 0 to K_) at positions k+1
												// e.g. alphabet size_ = 4 and K_ = 2: Y_ = 1 4 16 64
};

inline float** BackgroundModel::getV(){
    return v_;
}

#endif /* BACKGROUNDMODEL_H_ */
