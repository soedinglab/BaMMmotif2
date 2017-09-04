#ifndef GBaMM_H_
#define GBaMM_H_

#include <iostream>								// std::cout, std::endl
#include <vector>
#include <list>
#include <string>
#include <cstring> 								// std::strcpy
#include <numeric>								// std::iota
#include <iomanip>								// std::setprecision
#include <fstream>
#include <limits>								// std::numeric_limits<int>::max();
#include <random>
#include <algorithm>							// std::min

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


#include "../init/SequenceSet.h"
#include "../init/Alphabet.h"
#include "../getopt_pp/getopt_pp.h"                          // GetOpt function

class GBaMM{

public:

	static char*		outputDirectory; 		// output directory

	static char*		posSequenceFilename;	// filename of positive sequence FASTA file
	static std::string	posSequenceBasename;	// basename of positive sequence FASTA file
	static SequenceSet*	posSequenceSet;			// positive sequence set
	static bool 		maskPosSequenceSet;		// mask motif patterns from positive sequence set

	static char*		negSequenceFilename;	// filename of negative sequence FASTA file
	static std::string	negSequenceBasename;	// basename of negative sequence FASTA file
	static SequenceSet*	negSequenceSet;			// negative sequence set
	static bool			negSeqGiven;			// a flag for the negative sequence given by users

	// weighting options
	static char*		intensityFilename;		// filename of intensity file (i.e. for HT-SELEX data)

	// sequence set options
	static char* 		alphabetType;			// provide alphabet type
	static bool			ss;						// only search on single strand sequences

	// initial model(s) options
	static char*		initialModelFilename;	// filename of initial model
	static std::string	initialModelBasename;	// basename of initial model
	static std::string	initialModelTag;		// tag for initializing the model
	static size_t		num;					// number of init that are to be optimized
	static bool			mops;					// learn MOPS model
	static bool			zoops;					// learn ZOOPS model

	// model options
	static size_t		modelOrder;				// model order
	static std::vector<float> modelAlpha;		// initial alphas
	static float		modelBeta;				// alpha_k = beta x gamma^k for k > 0
	static float		modelGamma;
	static std::vector<size_t>	addColumns;		// add columns to the left and right of init used to initialize Markov init
    static bool			interpolate;				// calculate prior probabilities from lower-order probabilities
    											// instead of background frequencies of mono-nucleotides
    static bool			interpolateBG;			// calculate prior probabilities from lower-order probabilities
    											// instead of background frequencies of mono-nucleotides
	static float		q;						// prior probability for a positive sequence to contain a motif

	// background model options
    static char*		bgModelFilename;		// path to the background model file
    static bool			bgModelGiven;			// flag to show if the background model is given or not
	static size_t		bgModelOrder;			// background model order, defaults to 2
	static std::vector<float> bgModelAlpha;		// background model alpha

	// EM options
	static bool			EM;						// flag to trigger EM learning

	// Gibbs sampling options
	static bool			CGS;					// flag to trigger Collapsed Gibbs sampling
	static bool			noInitialZ;				// enable initializing z with one E-step
	static bool			noAlphaOptimization;	// disable alpha optimization in CGS
	static bool			GibbsMHalphas;			// enable alpha sampling in CGS using Gibbs Metropolis-Hastings
	static bool			dissampleAlphas;		// enable alpha sampling in CGS using discretely sampling
	static bool			noZSampling;			// disable sampling of z in CGS
	static bool			noQSampling;			// disable sampling of q in CGS
	static bool			debugAlphas;

	// FDR options
	static bool			FDR;					// triggers False-Discovery-Rate (FDR) estimation
	static size_t		mFold;					// number of negative sequences as multiple of positive sequences
	static size_t		cvFold;					// number of cross-validation (cv) folds
	static size_t		sOrder;					// k-mer order for sampling negative sequence set

	// motif occurrence options
	static bool         scoreSeqset;			// write logOdds Scores of positive sequence set to disk
	static float        scoreCutoff;			// cutoff for logOdds scores to print out as motif hits

	// other options
	static bool			verbose;				// verbose printouts, defaults to false
	static bool         debugMode;				// verbose printouts for debugging, defaults to false
	static bool			saveBaMMs;				// write optimized BaMM(s) to disk
	static bool			savePRs;				// write the precision, recall, TP and FP
	static bool			savePvalues;			// write p-values for each log odds score from sequence set
	static bool			saveLogOdds;			// write the log odds of positive and negative sets to disk
	static bool			saveInitialBaMMs;		// write out the initial model to disk
	static bool 		saveBgModel;			// write out background model to disk
	static bool			generatePseudoSet;		// test for alpha learning

	// flags for developers
	static bool			makeMovie;				// print out bamms in each iteration while optimizing
	static bool 		optimizeQ;				// optimize hyperparameter q in EM algorithm
	static void         debug();

	static void         init( int nargs, char* args[] );
	static void         destruct();
	static void			printStat();
	static char* 		String( const char *s );// convert const char* to string, for GetOpt library

	static std::mt19937	rngx;

private:

	static int	        readArguments( int nargs, char* args[] );
	static void	        printHelp();
};

// format cast for GetOpt function
inline char* GBaMM::String( const char *s ){
	return strdup( s );
}

namespace GetOpt{
	template <> inline _Option::Result convert<char*>( const std::string& s,
			char*& d, std::ios::fmtflags ){
		_Option::Result ret = _Option::BadType;
		d = GBaMM::String( s.c_str() );
		ret = _Option::OK;
		return ret;
	}
}


#endif /* GBaMM_H_ */
