/*
 * Global.cpp
 *
 *  Created on: Apr 4, 2016
 *      Author: wanwan
 */

#include <sys/stat.h> // get file status

#include "Global.h"

char* 				Global::outputDirectory = NULL;        	// output directory

char* 				Global::posSequenceFilename = NULL;		// filename of positive sequence fasta file
char* 				Global::negSequenceFilename = NULL;		// filename of negative sequence fasta file

char const*			Global::alphabetString = "ACGT";		// defaults to ACGT but may later be extended to ACGTH(hydroxymethylcytosine) or similar
bool				Global::revcomp = false;				// also search on reverse complement of sequences

SequenceSet* 		Global::posSequenceSet = NULL;			// positive Sequence Set
SequenceSet* 		Global::negSequenceSet = NULL;			// negative Sequence Set

char* 				Global::intensityFilename = NULL;		// filename of intensity file (i.e. for HT-SELEX data)
// further weighting options...

// files to initialize model(s)
char* 				Global::GIMMEpatternFilename = NULL;	// filename of GIMMEpattern file
char* 				Global::bindingSitesFilename = NULL;	// filename of binding sites file
char* 				Global::PWMFilename = NULL;				// filename of PWM file
char* 				Global::iIMMFilename = NULL;			// filename of Markov model (.iimm) file

// model options
unsigned int 		Global::modelOrder = 2;					// model order
float** 			Global::modelAlpha;						// initial alphas
std::vector<int>	Global::addColumns( 0, 0 );				// add columns to the left and right of models used to initialize Markov models
bool				Global::noLengthOptimization = false;	// disable length optimization

// background model options
unsigned int 		Global::bgModelOrder = 2;				// background model order, defaults to 2
float** 			Global::bgModelAlpha;					// background model alphas

// EM options
unsigned int		Global::maxEMIterations;				// maximum number of iterations
float 				Global::epsilon = 0.001f;				// likelihood convergence parameter

bool				Global::noAlphaOptimization = false;	// disable alpha optimization
bool				Global::noQOptimization = false;		// disable q optimization

// FDR options
bool				Global::FDR = false;					// triggers False-Discovery-Rate (FDR) estimation
unsigned int 		Global::mFold = 20;						// number of negative sequences as multiple of positive sequences
unsigned int 		Global::nFolds = 4;						// number of cross-validation folds
std::vector< std::vector<int> > Global::posFoldIndices;		// sequence indices for each cross-validation fold
std::vector< std::vector<int> > Global::negFoldIndices;		// sequence indices for each cross-validation fold
// further FDR options...

bool				Global::verbose = false;				// verbose printouts

int Global::readArguments( int nargs, char* args[] ){

	/*
	 * Process input arguments:
	 * 1. Essential parameters
	 * 		* output directory: the first argument
	 * 		* positive sequence file: the second argument
	 * 2. flags
	 * 		* initial motif model
	 * 			* GIMMEpattern file
	 * 			* binding sites file
	 * 			* PWM file
	 * 			* iIMM file
	 * 			( They must not be provided simultaneously! )
	 * 		* model and background model options
	 * 			* modelOrder / bgModelOrder
	 * 			* modelAlpha / bgModelAlpha
	 * 		* --fdr
	 * 			* mfold: for generating negative sequence file
	 * 			* nfolds: for cross-validation
	 * 		* alphabetString
	 * 			* if NULL, by default: ACGT
	 *
	 */

	if( nargs < 3 ) {		// At least 2 parameters are required:
		printHelp();		// 1.outputDirectory; 2.posSequenceFilename.
		exit( -1 );
	}

	outputDirectory = args[1];
	createDirectory( outputDirectory );

	posSequenceFilename = args[2];

	for( int i = 3; i < nargs; i++ ){

		if( strcmp( args[i], "--alphabetString" ) == 0 ){
			alphabetString = args[i+1];
		}else if( strcmp( args[i], "--revcomp" ) == 0 ){
			revcomp = true;
		}else if( strcmp( args[i], "--negSequenceFile" ) == 0 ){
			negSequenceFilename = args[i+1];
		}else if( strcmp( args[i], "--intensityFile" ) == 0 ){
			intensityFilename = args[i+1];
		}else if( strcmp( args[i], "--GIMMEpatternFile" ) == 0 ){
			GIMMEpatternFilename = args[i+1];
		}else if( strcmp( args[i], "--bindingSitesFile" ) == 0 ){
			bindingSitesFilename = args[i+1];
		}else if( strcmp( args[i], "--PWMFile" ) == 0 ){
			PWMFilename = args[i+1];
		}else if( strcmp( args[i], "--iIMMFile" ) == 0 ){
			iIMMFilename = args[i+1];
		}else if( strcmp( args[i], "--modelOrder" ) == 0 ){
			modelOrder = atoi( args[i+1] );
		}else if( strcmp( args[i], "--modelAlpha" ) == 0 ){
			// ...
		}else if( strcmp( args[i], "--addColumns" ) == 0 ){
			int left = atoi( args[i+1] ), right = atoi( args[i+2] );
			// addColumns( left, right );
		}else if( strcmp( args[i], "--noLengthOptimization" ) == 0 ){
			noLengthOptimization = true;
		}else if( strcmp( args[i], "--bgModelOrder" ) == 0 ){
			bgModelOrder = atoi( args[i+1] );
		}else if( strcmp( args[i], "--bgModelAlpha" ) == 0 ){
			// ...
		}else if( strcmp( args[i], "--maxEMIterations" ) == 0 ){
			maxEMIterations = atoi( args[i+1] );
		}else if( strcmp( args[i], "--epsilon" ) == 0 ){
			epsilon = atof( args[i+1] );
		}else if( strcmp( args[i], "--noAlphaOptimization" ) == 0 ){
			noAlphaOptimization = true;
		}else if( strcmp( args[i], "--noQOptimization" ) == 0 ){
			noQOptimization = true;
		}else if( strcmp( args[i], "--fdr" ) == 0 ){
			FDR = true;
			mFold = atoi( args[i+1] );
			nFolds = atoi( args[i+2] );
		}else if( strcmp( args[i], "--verbose" ) == 0 ){
			verbose = true;
		}
	}
	return 0;
};

void Global::printHelp(){

};

void Global::createDirectory( const char* os ){

	/*
	 * Create output directory
	 * works for Unix-like operating systems
	 */
	struct stat fileStatus;
	if( stat( os, &fileStatus ) != 0 ){
		fprintf( stderr, "Output directory does not exist. "
				"New directory is created automatically.\n" );
		char* command = (char*)calloc( 1024, sizeof(char) );
		sprintf( command, "mkdir %s", os );
		if( system( command ) != 0 ){
			fprintf( stderr, "Directory %s could not be created.\n", os );
			exit( -1 );
		}
		free( command );
	}

};

void Global::generateFolds(){
	// generate posFoldIndices
	if( negSequenceFilename == NULL ){
		// assign reference to posFoldIndices
	} else{
		// generate negFoldIndices
	}
};



