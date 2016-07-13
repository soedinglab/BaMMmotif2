/*
 * Global.cpp
 *
 *  Created on: Apr 4, 2016
 *      Author: wanwan
 */

#include <sys/stat.h>   			// get file status

#include "Global.h"

char*               Global::outputDirectory = NULL;			// output directory

char*               Global::posSequenceFilename = NULL;		// filename of positive sequence FASTA file
char*               Global::negSequenceFilename = NULL;		// filename of negative sequence FASTA file

char*				Global::posSequenceBasename = NULL;		// basename of positive sequence FASTA file
char*				Global::negSequenceBasename = NULL;		// basename of negative sequence FASTA file
char*				Global::alphabetType = NULL;			// alphabet type is defaulted to standard which is ACGT
bool                Global::revcomp = false;				// also search on reverse complement of sequences

SequenceSet*        Global::posSequenceSet = NULL;			// positive Sequence Set
SequenceSet*        Global::negSequenceSet = NULL;			// negative Sequence Set

char*               Global::intensityFilename = NULL;		// filename of intensity file (i.e. for HT-SELEX data)
// further weighting options...

// files to initialize model(s)
char*               Global::BaMMpatternFilename = NULL;		// filename of BaMMpattern file
char*               Global::bindingSitesFilename = NULL;	// filename of binding sites file
char*               Global::PWMFilename = NULL;				// filename of PWM file
char*               Global::BaMMFilename = NULL;			// filename of Markov model (.bmm) file
char*				Global::initialModelBasename = NULL;	// basename of initial model

// model options
int        			Global::modelOrder = 2;					// model order
std::vector<float> 	Global::modelAlpha;						// initial alphas
std::vector<int>    Global::addColumns(2);					// add columns to the left and right of models used to initialize Markov models

// background model options
int        			Global::bgModelOrder = 2;				// background model order, defaults to 2
float				Global::bgModelAlpha = 10.0f;			// background model alpha

// EM options
unsigned int        Global::maxEMIterations = std::numeric_limits<int>::max();	// maximum number of iterations
//unsigned int        Global::maxEMIterations = 5;			// maximum number of iterations
float               Global::epsilon = 0.001f;				// threshold for likelihood convergence parameter
bool                Global::noAlphaOptimization = false;	// disable alpha optimization
bool                Global::noQOptimization = false;		// disable q optimization

// FDR options
bool                Global::FDR = false;					// triggers False-Discovery-Rate (FDR) estimation
unsigned int        Global::mFold = 20;						// number of negative sequences as multiple of positive sequences
unsigned int        Global::cvFold = 5;						// size of cross-validation folds
std::vector< std::vector<int> > Global::posFoldIndices; 	// sequence indices for each cross-validation fold
std::vector< std::vector<int> > Global::negFoldIndices; 	// sequence indices for each cross-validation fold
// further FDR options...

bool                Global::verbose = false;            	// verbose printouts

int* 				Global::powA = NULL;

void Global::init( int nargs, char* args[] ){

	readArguments( nargs, args );

	Alphabet::init( alphabetType );

	// read in positive sequence set
	posSequenceSet = new SequenceSet( posSequenceFilename );

	// read in or generate negative sequence set
	if( negSequenceFilename == NULL )						// use positive for negative sequence set
//		negSequenceSet = new SequenceSet( posSequenceFilename );
		negSequenceSet = new SequenceSet( *posSequenceSet );
	else													// read in negative sequence set
		negSequenceSet = new SequenceSet( negSequenceFilename );

	if( verbose ){
		std::cout << "For positive set:	" << posSequenceSet->getN() << " sequences. "
				"max.length: " << posSequenceSet->getMaxL() << ", min.length: " <<
				posSequenceSet->getMinL() << std::endl << "			base frequencies: " ;
		for( int i = 0; i < Alphabet::getSize(); i++ )
			std::cout << std::scientific << Alphabet::getAlphabet()[i]
			          << " " << posSequenceSet->getBaseFrequencies()[i] << ", ";
		std::cout << std::endl;

		std::cout << "For negative set:	" << negSequenceSet->getN() << " sequences. "
				"max.length: " << negSequenceSet->getMaxL() << ", min.length: " <<
				negSequenceSet->getMinL() << std::endl << "			base frequencies: " ;
		for( int i = 0; i < Alphabet::getSize(); i++ )
			std::cout << std::scientific << Alphabet::getAlphabet()[i]
			          << " " << negSequenceSet->getBaseFrequencies()[i] << ", ";
		std::cout << std::endl;
	}

	// generate folds (fill posFoldIndices and negFoldIndices)
	generateFolds( posSequenceSet->getN(), negSequenceSet->getN(), cvFold );

	// optional: read in sequence intensities (header and intensity columns?)
	if( intensityFilename != 0 ){
		;// read in sequence intensity
	}

	powA = new int[modelOrder+2];
	for( int k = 0; k < modelOrder+2; k++ )
		powA[k] = ipow( Alphabet::getSize(), k );
}

int Global::readArguments( int nargs, char* args[] ){

	/*
	 * Process input arguments:
	 * 1. Essential parameters:
	 * 		* output directory: the first argument
	 * 		* positive sequence file: the second argument
	 * 2. flags
	 * 		* initial motif model
	 * 			* BaMMpattern file
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
	 * 		* alphabetType
	 * 			* if NULL, by default: ACGT
	 *
	 */

	if( nargs < 3 ) {			// At least 2 parameters are required, but this is not precise, needed to be specified!!
		fprintf( stderr, "Error: Arguments are missing! \n" );
		printHelp();
		exit( -1 );
	}

	// read in the output directory and create it
	outputDirectory = args[1];
	createDirectory( outputDirectory );

	// read in the positive sequence file
	posSequenceFilename = args[2];
	posSequenceBasename = baseName( posSequenceFilename );

	/**
	 * read command line to get options
	 * process flags from user
	 */

	GetOpt::GetOpt_pp opt( nargs, args );

	// negative sequence set
	if( opt >> GetOpt::OptionPresent( "negSeqFile" ) ){
		opt >> GetOpt::Option( "negSeqFile", negSequenceFilename );
		negSequenceBasename = baseName( negSequenceFilename );
	}

	// Alphabet Type
	if( opt >> GetOpt::OptionPresent( "alphabet" ) ){
		opt >> GetOpt::Option( "alphabet", alphabetType );
	} else {
		alphabetType = new char[9];
		strcpy( alphabetType, "STANDARD");
	}

	opt >> GetOpt::OptionPresent( "revcomp", revcomp );

	// for HT-SLEX data
	opt >> GetOpt::Option( "intensityFile", intensityFilename );

	// get initial motif files
	if( opt >> GetOpt::OptionPresent( "BaMMpatternFile" ) ){
		opt >> GetOpt::Option( "BaMMpatternFile", BaMMpatternFilename );
		initialModelBasename = baseName( BaMMpatternFilename );
	} else if ( opt >> GetOpt::OptionPresent( "bindingSitesFile" ) ){
		opt >> GetOpt::Option( "bindingSitesFile", bindingSitesFilename );
		initialModelBasename = baseName( bindingSitesFilename );
	} else if ( opt >> GetOpt::OptionPresent( "PWMFile" ) ){
		opt >> GetOpt::Option( "PWMFile", PWMFilename );
		initialModelBasename = baseName( PWMFilename );
	} else if( opt >> GetOpt::OptionPresent( "BMMFile" ) ){
		opt >> GetOpt::Option( "BMMFile", BaMMFilename );
		initialModelBasename = baseName( BaMMFilename );
	} else {
		fprintf( stderr, "Error: No initial model is provided.\n" );
		exit( -1 );
	}

	// model options
	opt >> GetOpt::Option( 'k', "modelOrder", modelOrder );
	if( opt >> GetOpt::OptionPresent( "modelAlpha" ) || opt >> GetOpt::OptionPresent( 'a' ) )
		opt >> GetOpt::Option( 'a', "modelAlpha", modelAlpha );
	else 		// defaults modelAlpha.at(k) to 10.0
		for( int k = 0; k < modelOrder+1; k++ )
			modelAlpha.push_back( 10.0f );

	if( opt >> GetOpt::OptionPresent( "addColumns" ) ){
		addColumns.clear();
		opt >> GetOpt::Option( "addColumns", addColumns );
		if( addColumns.size() < 1 || addColumns.size() > 2 ){
			fprintf( stderr, "--addColumns format error.\n" );
			exit( -1 );
		}
		if( addColumns.size() == 1 )
			addColumns.resize( 2, addColumns.back() );
	} else {
		addColumns.at(0) = 0;
		addColumns.at(1) = 0;
	}

	// background model options
	opt >> GetOpt::Option( 'K', "bgModelOrder", bgModelOrder );
	opt >> GetOpt::Option( 'A', "bgModelAlpha", bgModelAlpha );

	// em options
	opt >> GetOpt::Option( "maxEMIterations",  maxEMIterations );
	opt >> GetOpt::Option( "epsilon", epsilon );
	opt >> GetOpt::OptionPresent( "noAlphaOptimization", noAlphaOptimization );
	opt >> GetOpt::OptionPresent( "noQOptimization", noQOptimization );

	// FDR options
	if( opt >> GetOpt::OptionPresent( "FDR", FDR) ){
		opt >> GetOpt::Option( "m", mFold  );
		opt >> GetOpt::Option( "n", cvFold );
	}

	opt >> GetOpt::OptionPresent( "verbose", verbose );

	return 0;
}

void Global::createDirectory( char* dir ){
	struct stat fileStatus;
	if( stat( dir, &fileStatus ) != 0 ){
		fprintf( stderr, "Status: Output directory does not exist. "
				"New directory is created automatically.\n" );
		char* command = ( char* )calloc( 1024, sizeof( char ) );
		sprintf( command, "mkdir %s", dir );
		if( system( command ) != 0 ){
			fprintf( stderr, "Directory %s could not be created.\n", dir );
			exit( -1 );
		}
		free( command );
	}
}

char* Global::baseName( char* filepath ){
	std::string path = filepath;
	char* base;
	int i = 0, start = 0, end = 0;
	while( path[++i] != '\0' )
		if( path[i] == '.' )
			end = i - 1;
	while( --i != 0 && path[i] != '/' )
		;
	if( i == 0 )
		start = 0;
	else
		start = i + 1;
	base = ( char* )malloc( ( end-start+2 ) * sizeof( char ) );
	for( i = start; i <= end; i++ )
		base[i-start] = path[i];
	base[i-start] = '\0';
	return base;
}

void Global::generateFolds( unsigned int posN, unsigned int negN, unsigned int fold ){
	posFoldIndices.resize( fold );
	negFoldIndices.resize( fold );
	// generate posFoldIndices
	unsigned int i = 0;
	while( i < posN ){
		for( unsigned int j = 0; j < fold; j++ ){
			posFoldIndices[j].push_back( i );
			i++;
		}
	}
	// generate negFoldIndices
	if( negSequenceFilename == NULL ){
		negFoldIndices = posFoldIndices;
	} else {
		// generate negFoldIndices
		i = 0;
		while( i < negN ){
			for( unsigned int j = 0; j < fold; j++ ){
				negFoldIndices[j].push_back( i );
				i++;
			}
		}
	}
}

// power function for integers
int Global::ipow( unsigned int base, int exp ){
    int result = 1;
    while( exp ){
        if( exp & 1 )
        	result *= base;
        exp >>= 1;
        base *= base;
    }
    return result;
}

void Global::printHelp(){
	printf("\n=================================================================\n");
	printf("\n Usage: ./BaMMmotif OUTDIR SEQFILE [options] \n\n");
	printf("\t OUTDIR:  output directory for all results. \n");
	printf("\t SEQFILE: file with sequences from positive set in FASTA format. \n");
	printf("\n Options: \n");
	// To be fulfilled ...
	printf(" \n");
	printf(" \n");
	printf("\n=================================================================\n");
}

void Global::destruct(){
	Alphabet::destruct();
	if( powA ) delete[] powA;
	if( alphabetType ) delete[] alphabetType;
	if( posSequenceBasename ) free( posSequenceBasename );
	if( negSequenceBasename ) free( negSequenceBasename );
	if( initialModelBasename ) free( initialModelBasename );
	if( posSequenceSet ) delete posSequenceSet;
	if( negSequenceFilename ) delete negSequenceSet;
	std::cout << "Destruct() for Global class works fine. \n";
}
