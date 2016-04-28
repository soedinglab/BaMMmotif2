/*
 * Global.h
 *
 *  Created on: Apr 1, 2016
 *      Author: wanwan
 */

#ifndef GLOBAL_H_
#define GLOBAL_H_

#include <vector>
#include <string>
#include <iostream>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "SequenceSet.h"
#include "Alphabet.h"
#include "lib/getopt_pp.h"		// GetOpt function

class Global {

public:

	static char*		outputDirectory;      	// output directory

	static char*		posSequenceFilename;	// filename of positive sequence FASTA file
	static char*		negSequenceFilename;	// filename of negative sequence FASTA file

	static char*		posSequenceBasename;	// basename of positive sequence FASTA file
	static char*		negSequenceBasename;	// basename of negative sequence FASTA file


	static char* 		alphabetType;			// provide alphabet type
	static bool			revcomp;				// also search on reverse complement of sequences, defaults to false

	static SequenceSet*	posSequenceSet;			// positive Sequence Set
	static SequenceSet*	negSequenceSet;			// negative Sequence Set

	static char*		intensityFilename;		// filename of intensity file (i.e. for HT-SELEX data)
	// further weighting options...

	// files to initialize model(s)
	static char*		BaMMpatternFilename;	// filename of BaMMpattern file
	static char*		bindingSitesFilename;	// filename of binding sites file
	static char*		PWMFilename;			// filename of PWM file
	static char*		BMMFilename;			// filename of Markov model (.bmm) file

	// model options
	static unsigned int	modelOrder;				// model order
	static std::vector<float>	modelAlpha;		// initial alphas
	static std::vector<int>		addColumns;		// add columns to the left and right of models used to initialize Markov models
	static bool			noLengthOptimization;	// disable length optimization

	// background model options
	static unsigned int	bgModelOrder;			// background model order, defaults to 2
	static float		bgModelAlpha;	// background model alpha

	// EM options
	static unsigned int	maxEMIterations;		// maximum number of iterations
	static float		epsilon;				// likelihood convergence parameter

	static bool			noAlphaOptimization;	// disable alpha optimization
	static bool			noQOptimization;		// disable q optimization

	// FDR options
	static bool			FDR;					// triggers False-Discovery-Rate (FDR) estimation
	static unsigned int	mFold;					// number of negative sequences as multiple of positive sequences
	static unsigned int	nFolds;					// number of cross-validation folds
	static std::vector< std::vector<int> >	posFoldIndices;	// sequence indices for each cross-validation fold
	static std::vector< std::vector<int> >	negFoldIndices;	// sequence indices for each cross-validation fold
	// further FDR options...

	static bool			verbose;				// verbose printouts, defaults to false

	static void         init( int nargs, char* args[] );

	static void         destruct();

private:

	static int	        readArguments( int nargs, char* args[] );
	static void	        createDirectory( char* dir );
	static void	        generateFolds( unsigned int posN,unsigned int negN, unsigned int fold );
	static void	        printHelp();
};

#endif /* GLOBAL_H_ */
