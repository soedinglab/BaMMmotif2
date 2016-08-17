#ifndef BACKGROUNDMODEL_H_
#define BACKGROUNDMODEL_H_

#include <fstream>	// e.g. std::ofstream
#include <iomanip>	// e.g. std::setprecision
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
			         std::vector< std::vector<int> > foldIndices = std::vector< std::vector<int> >(),
			         std::vector<int> folds = std::vector<int>() );

	BackgroundModel( std::string filePath );

	~BackgroundModel();

	std::string getName();
	float**		getVbg();

	void 	print();
	void 	write( char* dir );

private:

	void 	calculateVbg();  	// calculate conditional probabilities from counts

	std::string			name_;	// basename of sequence set file
	int**				n_bg_;	// oligomer counts
	float** 			v_bg_;	// oligomer conditional probabilities
	int					K_;		// order
	std::vector<float>	A_;		// order-specific alphas

	std::vector<int>	Y_;		// contains 1 at position 0
								// and the number of oligomers y for increasing order k (from 0 to K_) at positions k+1
								// e.g. alphabet size_ = 4 and K_ = 2: Y_ = 1 4 16 64
};

#endif /* BACKGROUNDMODEL_H_ */
