/*
 * Global.h
 *
 *  Created on: Apr 1, 2016
 *      Author: wanwan
 */

#ifndef GLOBAL_H_
#define GLOBAL_H_

#include "SequenceSet.h"

class Global {

public:

	static char* 			outputDirectory;        	// output directory

	static char*			alphabetString;				// defaults to ACGT but may later be extended to ACGTH(ydroxymethylcytosine) or similar
	static char* 			posSequenceFilename;		// filename of positive sequence fasta file
	static char* 			negSequenceFilename;		// filename of negative sequence fasta file

	static bool				revcomp;					// also search on reverse complement of sequences

	static SequenceSet* 	posSequenceSet;				// positive SequenceSet
	static SequenceSet* 	negSequenceSet;				// negative SequenceSet

	static char* 			intensityFilename;			// filename of intensity file (i.e. for HT-SELEX data)
	// further weighting options...

	// files to initialize model(s)
	static char* 			GIMMEpatternFilename;		// filename of GIMMEpattern file
	static char* 			bindingSitesFilename;		// filename of binding sites file
	static char* 			PWMFilename;				// filename of PWM file
	static char* 			iIMMFilename;				// filename of Markov model (.iimm) file

	// model options
	static int 				modelOrder;					// model order
	static float** 			modelAlpha;					// initial alphas
	static std::vector<int> addColumns;					// add columns to the left and right of models used to initialize Markov models
	static bool				noLengthOptimization;		// disable length optimization

	// background model options
	static int 				bgModelOrder;				// background model order
	static float** 			bgModelAlpha;				// background model alphas

	// EM options
	static int				maxEMIterations;			// maximum number of iterations
	static float 			epsilon;					// likelihood convergence parameter

	static bool				noAlphaOptimization;		// disable alpha optimization
	static bool				noQOptimization;			// disable q optimization

	// FDR options
	static int 				M;							// number of negative sequences as multiple of positive sequences
	static int 				nFolds;						// number of cross-validation folds
	static std::vector<std::vector<int>> posFoldIndices;// sequence indices for each cross-validation fold
	static std::vector<std::vector<int>> negFoldIndices;// sequence indices for each cross-validation fold

	// further FDR options...

	static bool				verbose;					// verbose printouts

	static void init( int nargs, char* args[] ){
		readArguments( nargs, args );
		Alphabet::init( alphabetString );
		// read in positive (and negative) sequence set
		// generate folds (fill posFoldIndices and negFoldIndices)
		// optional: read in sequence intensities (header and intensity columns?)
	}

	static void destruct(){
		Alphabet::destruct();
		// ...
	}

private:

	static int  readArguments( int nargs, char *args[] );
	static void printHelp();

	static void createDirectory();
	static void generateFolds(){
		// generate posFoldIndices
		if( negSequenceFilename == NULL ){
			// assign reference to posFoldIndices
		} else{
			// generate negFoldIndices
		}
	};
};

#endif /* GLOBAL_H_ */
