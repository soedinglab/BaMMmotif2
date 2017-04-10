#ifndef GLOBAL_H_
#define GLOBAL_H_

#include <iostream>											// std::cout, std::endl
#include <vector>
#include <list>
#include <string>
#include <cstring> 											// std::strcpy
#include <numeric>											// std::iota
#include <iomanip>											// std::setprecision
#include <fstream>
#include <limits>											// std::numeric_limits<int>::max();
#include <random>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <algorithm>                                        // std::min

#include "../shared/SequenceSet.h"
#include "../shared/Alphabet.h"
#include "../getopt_pp/getopt_pp.h"							// GetOpt function

class Global{

public:

	static char*		outputDirectory; 					// output directory

	static char*		posSequenceFilename;				// filename of positive sequence FASTA file
	static std::string	posSequenceBasename;				// basename of positive sequence FASTA file
	static SequenceSet*	posSequenceSet;						// positive Sequence Set
	static std::vector<std::vector<int>> posFoldIndices;	// sequence indices for positive sequence set

	static char*		negSequenceFilename;				// filename of negative sequence FASTA file
	static std::string	negSequenceBasename;				// basename of negative sequence FASTA file
	static SequenceSet*	negSequenceSet;						// negative Sequence Set
	static std::vector<std::vector<int>> negFoldIndices;	// sequence indices for given negative sequence set
	static bool			negSeqGiven;						// a flag for the negative sequence given by users

	// weighting options
	static char*		intensityFilename;					// filename of intensity file (i.e. for HT-SELEX data)

	// sequence set options
	static char* 		alphabetType;						// provide alphabet type
	static bool			ss;									// only search on single strand sequences

	// initial model(s) options
	static char*		BaMMpatternFilename;				// filename of BaMMpattern file
	static char*		bindingSiteFilename;				// filename of binding sites file
	static char*		PWMFilename;						// filename of PWM file
	static char*		BaMMFilename;						// filename of Markov model (.bamm) file
	static std::string	initialModelBasename;				// basename of initial model
	static int			num;								// number of models that are to be optimized
	static bool			mops;								// learn MOPS model
	static bool			zoops;								// learn ZOOPS model

	// model options
	static int			modelOrder;							// model order
	static std::vector<float> modelAlpha;					// initial alphas
	static float		modelBeta;							// alpha_k = beta x gamma^k for k > 0
	static float		modelGamma;
	static std::vector<int>	addColumns;						// add columns to the left and right of models used to initialize Markov models
    static bool        interpolate;                         // calculate prior probabilities from lower-order probabilities
                                                            // instead of background frequencies of mono-nucleotides
    static bool        interpolateBG;                       // calculate prior probabilities from lower-order probabilities
                                                            // instead of background frequencies of mono-nucleotides

	// background model options
	static int			bgModelOrder;						// background model order, defaults to 2
	static std::vector<float> bgModelAlpha;					// background model alpha

	// EM options
	static bool			EM;									// flag to trigger EM learning
	static int			maxEMIterations;					// maximum number of iterations for EM
	static float		epsilon;							// threshold for likelihood convergence parameter
	static bool			noQOptimization;					// disable q optimization
	static float		q;									// prior probability for a positive sequence to contain a motif

	// CGS (Collapsed Gibbs sampling) options
	static bool			CGS;								// flag to trigger Collapsed Gibbs sampling
	static int 			maxCGSIterations;					// maximum number of iterations for CGS
	static bool			noInitialZ;							// enable initializing z with one E-step
	static bool			noAlphaOptimization;				// disable alpha optimization in CGS
	static bool			alphaSampling;						// enable alpha sampling in CGS
	static bool			noQSampling;						// disable q sampling in CGS
	static float		eta;								// learning rate for optimizing alphas, only for tuning
	static int			interval;							// interval for sampling z and q, only for tuning
	static bool			debugAlphas;

	// FDR options
	static bool			FDR;								// triggers False-Discovery-Rate (FDR) estimation
	static int			mFold;								// number of negative sequences as multiple of positive sequences
	static int			cvFold;								// number of cross-validation (cv) folds
	static int 			sOrder;								// the k-mer order for sampling negative sequence set

	// scoring options
	static float        scoreCutoff;                        // cutoff for logOdds scores to print out as motif hits

	// other options
	static bool			verbose;							// verbose printouts, defaults to false
	static bool         debugMode;                          // verbose printouts for debugging, defaults to false
	static bool			saveBaMMs;							// write optimized BaMM(s) to disk
	static bool			savePRs;							// write the precision, recall, TP and FP
	static bool			savePvalues;						// write p-values for each log odds score from sequence set
	static bool			saveLogOdds;						// write the log odds of positive and negative sets to disk
	static bool			saveInitialModel;					// write out the initial model to disk
	static bool 		saveBgModel;						// write out background model to disk
	static bool         scoreSeqset;                        // write logOdds Scores of positive sequence set to disk
	static int			Yk;									// the counts of numbers in Y_ array
	static bool			generatePseudoSet;					// test for alpha learning

	static void         init( int nargs, char* args[] );
	static void         destruct();
	static void         debug();

	static char* 		String( const char *s );			// convert const char* to string, for GetOpt library

	static std::mt19937	rngx;

private:

	static int	        readArguments( int nargs, char* args[] );
	static void	        printHelp();
};

// format cast for GetOpt function
inline char* Global::String( const char *s ){
	return strdup( s );
}

namespace GetOpt{
	template <> inline _Option::Result convert<char*>( const std::string& s, char*& d, std::ios::fmtflags ){
		_Option::Result ret = _Option::BadType;
		d = Global::String( s.c_str() );
		ret = _Option::OK;
		return ret;
	}
}


#endif /* GLOBAL_H_ */
