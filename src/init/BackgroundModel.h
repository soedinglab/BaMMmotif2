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
#include "../Global/utils.h"

class BackgroundModel{

public:

	BackgroundModel( std::vector<Sequence*> sequenceSet,
					size_t order,
			        std::vector<float> alpha,
			        bool interpolate = true,
					std::string basename = "" );

	BackgroundModel( std::string filePath );
    BackgroundModel( char* filePath , int K, float A );

    ~BackgroundModel();

	std::string getName();
	size_t		getOrder();
	float**		getV();

	void 		expV();
	void 		logV();

	bool		vIsLog();

	// calculate the log likelihood for the sequence set
	// afterwards v_ contains log probabilities
	double 	calculateLogLikelihood( std::vector<Sequence*> sequenceSet );

	// calculate positional likelihoods for the sequence set
	// and write likelihoods to file
	void 	calculatePosLikelihoods( std::vector<Sequence*> sequenceSet,
			                      	  char* odir  );


	void 	print();
    // afterwards v_ contains exp probabilities
	void 	write( char* dir, std::string basename );

private:

	// calculate probabilities from conditional probabilities
	void	calculateProbabilities( float** p );
	// calculate conditional probabilities from counts
	void 	calculateV();

	std::string			basename_;			// basename of sequence set file

	size_t**			n_;					// oligomer counts
	float** 			v_;					// oligomer conditional probabilities

	bool				vIsLog_ = false;	// v_ contains log probabilities

	size_t				K_;					// order
	std::vector<float>	A_;					// order-specific alphas

	bool				interpolate_ = true;// calculate prior probabilities
											// from lower-order probabilities
											// instead of background frequencies
											// of mononucleotides

	std::vector<size_t>	Y_;
};

inline float** BackgroundModel::getV(){
    return v_;
}

#endif /* BACKGROUNDMODEL_H_ */
