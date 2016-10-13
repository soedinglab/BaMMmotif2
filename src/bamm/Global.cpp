#include <sys/stat.h>   			// get file status

#include "Global.h"

char*               Global::outputDirectory = NULL;			// output directory

char*               Global::posSequenceFilename = NULL;		// filename of positive sequence FASTA file
char*				Global::posSequenceBasename = NULL;		// basename of positive sequence FASTA file
SequenceSet*        Global::posSequenceSet = NULL;			// positive Sequence Set
std::vector< std::vector<int> >	Global::posFoldIndices;		// sequence indices for positive sequence set

char*               Global::negSequenceFilename = NULL;		// filename of negative sequence FASTA file
char*				Global::negSequenceBasename = NULL;		// basename of negative sequence FASTA file
SequenceSet*        Global::negSequenceSet = NULL;			// negative Sequence Set
std::vector< std::vector<int> >	Global::negFoldIndices;		// sequence indices for given negative sequence set
bool				Global::negGiven = false;				// a flag for the given negative sequence set by user

// weighting options
char*               Global::intensityFilename = NULL;		// filename of intensity file (i.e. for HT-SELEX data)

char*				Global::alphabetType = NULL;			// alphabet type is defaulted to standard which is ACGT
bool                Global::revcomp = false;				// also search on reverse complement of sequences

// initial model(s) options
char*               Global::BaMMpatternFilename = NULL;		// filename of BaMMpattern file
char*               Global::bindingSiteFilename = NULL;		// filename of binding sites file
char*               Global::PWMFilename = NULL;				// filename of PWM file
char*               Global::BaMMFilename = NULL;			// filename of Markov model (.bamm) file
char*				Global::initialModelBasename = NULL;	// basename of initial model

// model options
int        			Global::modelOrder = 2;					// model order
std::vector<float> 	Global::modelAlpha( modelOrder+1, 1.0f );	// initial alphas
float				Global::modelBeta = 20.0f;				// alpha_k = beta x gamma^(k-1) for k > 0
float				Global::modelGamma = 3.0f;
std::vector<int>    Global::addColumns(2);					// add columns to the left and right of initial model
bool                Global::interpolate = true;             // calculate prior probabilities from lower-order probabilities
                                                            // instead of background frequencies of mononucleotides
bool                Global::interpolateBG = true;           // calculate prior probabilities from lower-order probabilities
                                                            // instead of background frequencies of mononucleotides

// background model options
int        			Global::bgModelOrder = 2;				// background model order, defaults to 2
std::vector<float>	Global::bgModelAlpha( bgModelOrder+1, 1.0f );	// background model alpha

// EM options
unsigned int        Global::maxEMIterations = std::numeric_limits<int>::max();	// maximum number of iterations
float               Global::epsilon = 0.001f;				// threshold for likelihood convergence parameter
bool                Global::noAlphaOptimization = false;	// disable alpha optimization
bool                Global::noQOptimization = false;		// disable q optimization
bool				Global::setSlow = false;				// develop with the slow EM version
bool				Global::logEM = false;					// calculation EM steps in log space

// FDR options
bool                Global::FDR = false;					// triggers False-Discovery-Rate (FDR) estimation
int        			Global::mFold = 20;						// number of negative sequences as multiple of positive sequences
int        			Global::cvFold = 5;						// size of cross-validation folds
// further FDR options...

// printout options
bool                Global::verbose = false;
bool                Global::debugMode = false;              // debug-mode: prints out everything.
bool                Global::saveInitBaMMs = false;
bool				Global::saveBaMMs = true;

void Global::init( int nargs, char* args[] ){

	readArguments( nargs, args );

	Alphabet::init( alphabetType );

	// read in positive and negative sequence set
	posSequenceSet = new SequenceSet( posSequenceFilename );
	negSequenceSet = new SequenceSet( negSequenceFilename );

	// generate fold indices for positive and negative sequence set
	Global::posFoldIndices = generateFoldIndices( posSequenceSet->getN(), cvFold );
	Global::negFoldIndices = generateFoldIndices( negSequenceSet->getN(), cvFold );

	// optional: read in sequence intensities (header and intensity columns?)
	if( intensityFilename != 0 ){
		;// read in sequence intensity
	}
}

int Global::readArguments( int nargs, char* args[] ){

	/**
	 * read command line to get options
	 * process flags from user
	 */

	if( nargs < 3 ) {										// At least 2 parameters are required, TODO:but this is not precise, needed to be specified!!
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

	GetOpt::GetOpt_pp opt( nargs, args );

	if( opt >> GetOpt::OptionPresent( 'h', "help" ) ){
		printHelp();
		exit( -1 );
	}

	// negative sequence set
	if( opt >> GetOpt::OptionPresent( "negSeqFile" ) ){
		opt >> GetOpt::Option( "negSeqFile", negSequenceFilename );
	}else{
	    negSequenceFilename = posSequenceFilename;
	}
	negSequenceBasename = baseName( negSequenceFilename );

	// Alphabet Type
	if( opt >> GetOpt::OptionPresent( "alphabet" ) ){
		opt >> GetOpt::Option( "alphabet", alphabetType );
	} else {
		alphabetType = new char[9];
		strcpy( alphabetType, "STANDARD");
	}

	opt >> GetOpt::OptionPresent( "reverseComp", revcomp );

	// for HT-SELEX data
	opt >> GetOpt::Option( "intensityFile", intensityFilename );

	// get initial model files
	if( opt >> GetOpt::OptionPresent( "BaMMpatternFile" ) ){
		opt >> GetOpt::Option( "BaMMpatternFile", BaMMpatternFilename );
		initialModelBasename = baseName( BaMMpatternFilename );
	} else if ( opt >> GetOpt::OptionPresent( "bindingSiteFile" ) ){
		opt >> GetOpt::Option( "bindingSiteFile", bindingSiteFilename );
		initialModelBasename = baseName( bindingSiteFilename );
	} else if ( opt >> GetOpt::OptionPresent( "PWMFile" ) ){
		opt >> GetOpt::Option( "PWMFile", PWMFilename );
		initialModelBasename = baseName( PWMFilename );
	} else if( opt >> GetOpt::OptionPresent( "BMMFile" ) ){
		opt >> GetOpt::Option( "BaMMFile", BaMMFilename );
		initialModelBasename = baseName( BaMMFilename );
	} else {
		fprintf( stderr, "Error: No initial model is provided.\n" );
		exit( -1 );
	}

	// model options
	opt >> GetOpt::Option( 'k', "order", modelOrder );

	if( opt >> GetOpt::OptionPresent( 'a', "alpha" ) ){
		modelAlpha.clear();
		opt >> GetOpt::Option( 'a', "alpha", modelAlpha );
		if( static_cast<int>( modelAlpha.size() ) != modelOrder+1 ){
			if( static_cast<int>( modelAlpha.size() ) > modelOrder+1 ){
				modelAlpha.resize( modelOrder+1 );
			} else{
				modelAlpha.resize( modelOrder+1, modelAlpha.back() );
			}
		}
	} else {
		if( static_cast<int>( modelAlpha.size() ) != modelOrder+1 ){
			if( static_cast<int>( modelAlpha.size() ) > modelOrder+1 ){
				modelAlpha.resize( modelOrder+1 );
			} else{
				modelAlpha.resize( modelOrder+1, modelAlpha.back() );
			}
		}
		if( modelOrder > 0 ){
			for( int k = 1; k < modelOrder + 1; k++ ){
				modelAlpha[k] = modelBeta * powf( modelGamma, static_cast<float>( k ) - 1.0f );
			}
		}
	}

	if( opt >> GetOpt::OptionPresent( "extend" ) ){
		addColumns.clear();
		opt >> GetOpt::Option( "extend", addColumns );
		if( addColumns.size() < 1 || addColumns.size() > 2 ){
			fprintf( stderr, "--extend format error.\n" );
			exit( -1 );
		}
		if( addColumns.size() == 1 )
			addColumns.resize( 2, addColumns.back() );
	} else {
		addColumns.at(0) = 0;
		addColumns.at(1) = 0;
	}

	// background model options
	opt >> GetOpt::Option( 'K', "Order", bgModelOrder );
//	if( bgModelOrder > modelOrder ){
//		std::cerr << "The order of background model should not exceed the order of motif model!\n";
//		exit( -1 );
//	}
	if( opt >> GetOpt::OptionPresent( 'A', "Alpha" ) ){
		bgModelAlpha.clear();
		opt >> GetOpt::Option( 'A', "Alpha", bgModelAlpha );
		if( static_cast<int>( bgModelAlpha.size() ) != bgModelOrder+1 ){
			if( static_cast<int>( bgModelAlpha.size() ) > bgModelOrder+1 ){
				bgModelAlpha.resize( bgModelOrder+1 );
			} else{
				bgModelAlpha.resize( bgModelOrder+1, bgModelAlpha.back() );
			}
		}
	} else {
		if( static_cast<int>( bgModelAlpha.size() ) != bgModelOrder+1 ){
			if( static_cast<int>( bgModelAlpha.size() ) > bgModelOrder+1 ){
				bgModelAlpha.resize( bgModelOrder+1 );
			} else{
				bgModelAlpha.resize( bgModelOrder+1, bgModelAlpha.back() );
			}
		}
		if( bgModelOrder > 0 ){
			for( int k = 1; k < bgModelOrder + 1; k++ ){
				bgModelAlpha[k] = 10.0f;
//				bgModelAlpha[k] = modelBeta * powf( modelGamma, static_cast<float>( k ) - 1.0f );
			}
		}
	}
	// em options
	opt >> GetOpt::Option( "maxEMIterations", maxEMIterations );
	opt >> GetOpt::Option( 'e', "epsilon", epsilon );
	opt >> GetOpt::OptionPresent( "noAlphaOptimization", noAlphaOptimization );
	opt >> GetOpt::OptionPresent( "noQOptimization", noQOptimization );
	opt >> GetOpt::OptionPresent( "setSlow", setSlow );
	opt >> GetOpt::OptionPresent( "logEM", logEM );

	// FDR options
	if( opt >> GetOpt::OptionPresent( "FDR", FDR ) ){
		opt >> GetOpt::OptionPresent( "FDR", FDR );
		opt >> GetOpt::Option( 'm', "mFold", mFold  );
		opt >> GetOpt::Option( 'n', "cvFold", cvFold );
	}

	// printout options
	opt >> GetOpt::OptionPresent( "verbose", verbose );
	opt >> GetOpt::OptionPresent( "debug", debugMode );
	opt >> GetOpt::OptionPresent( "saveInitBaMMs", saveInitBaMMs );
	opt >> GetOpt::OptionPresent( "saveBaMMs", saveBaMMs );

	return 0;
}

void Global::printHelp(){
	printf("\n===============================================================================\n");
	printf("\n SYNOPSIS:	BaMMmotif OUTDIR SEQFILE [options] \n\n");
	printf("\t DESCRIPTION \n");
	printf("\t 		Learn Bayesian inhomogeneous Markov models (BaMM) from sequence Data.\n"
			"		The default extension of sequence file is .fasta\n\n");
	printf("\t OUTDIR:  output directory for all results. \n");
	printf("\t SEQFILE: file with sequences from positive set in FASTA format. \n\n");
	printf("\n OPTIONS: \n");
	printf("\n		Options for reading in sequence file: \n");
	printf("\n			--alphabet <STRING> \n"
			"				STANDARD.		For alphabet type ACGT; \n"
			"				METHYLC. 		For alphabet type ACGTM; \n"
			"				HYDROXYMETHYLC.	For alphabet type ACGTH; \n"
			"				EXTENDED.		For alphabet type ACGTMH. \n\n");
	printf("\n			--reverseComp \n"
			"				search motif on the reverse complementary sequence as well. \n\n");
	printf("\n		Options for HT-SELEX data: \n");
	printf("\n			--intensityFile	<STRING> \n"
			"				intensity file name. \n\n");
	printf("\n		Options for initial model: \n");
	printf("\n 			--BaMMpatternFile <STRING> \n"
			"				file that contains patterns. \n\n");
	printf("\n 			--bindingSitesFile <STRING> \n"
			"				file that contains binding sites. \n\n");
	printf("\n 			--PWMFile <STRING> \n"
			"				file that contains PWM data. \n");
	printf("\n 			--BaMMFile <STRING> \n"
			"				file that contains BaMM data. \n\n");
	printf("\n 		Options for inhomogeneous BaMM: \n");
	printf("\n 			-k, --order <INTEGER> \n"
			"				model Order. The default is 2. \n\n");
	printf("\n 			-a, --alpha <FLOAT> [<FLOAT>...] \n"
			"				Order-specific prior strength. The default is 1.0 (for k = 0) and\n"
			"				20 x 3^(k-1) (for k > 0). The options -b and -g are ignored.\n\n");
	printf("\n 			--extend <INTEGER>{1, 2} \n"
			"				Extend BaMMs by adding uniformly initialized positions to the left\n"
			"				and/or right of initial BaMMs. e.g. invoking with --extend 0 2 adds\n"
			"				two positions to the right of initial BaMMs. Invoking with --extend 2\n"
			"				adds two positions to both sides of initial BaMMs. By default, BaMMs\n"
			"				are not being extended.\n\n");
	printf("\n 		Options for homogeneous (background) BaMM: \n");
	printf("\n 			-K, --Order <INTEGER> \n"
			"				Order. The default is 2.\n"
			"				Order of background model should not exceed order of motif model.\n\n");
	printf("\n 			-A, --Alpha <FLOAT> \n"
			"				Prior strength. The default is 10.0.\n\n");
	printf("\n 		Options for EM: \n");
	printf("\n 			-q <FLOAT> \n"
			"				Prior probability for a positive sequence to contain a motif. The"
			"				default is 0.9.\n\n");
	printf("\n 			-e, --epsilon <FLOAT> \n"
			"				The EM algorithm is deemed to be converged when the sum over the\n"
			"				absolute differences in BaMM probabilities from successive EM rounds\n"
			"				is smaller than epsilon. The default is 0.001.\n\n");
	printf("\n 			--maxEMIterations <INTEGER> (*) \n"
			"				Limit the number of EM iterations. *For developers.\n\n");
	printf("\n 			--noAlphaOptimization (*) \n"
			"				disable alpha optimization. Defaults to false. *For developers.\n\n");
	printf("\n 			--noQOptimization (*) \n"
			"				disable q optimization. Defaults to false. *For developers.\n\n");
	printf("\n 			--setSlow (*)\n"
			"				use slow version of EM. Defaults to false. *For developers.\n\n");
	printf("\n 			--logEM (*)\n"
			"				calculate EM in log space. Defaults to false. *For developers.\n\n");
	printf("\n 		Options for FDR: \n");
	printf("\n 			--FDR \n"
			"				triggers False-Discovery-Rate (FDR) estimation.\n\n");
	printf("\n 			-m, --mFold <INTEGER>\n"
			"				number of negative sequences as multiple of positive sequences."
			"				The default is 20.\n\n");
	printf("\n 			-n, --cvFold <INTEGER>\n"
			"				number of cross-validation folds. The default is 5.\n\n");
	printf("\n 		Options for output:	\n");
	printf("\n 			--verbose \n"
			"				verbose printouts. Defaults to false.\n\n");
	printf("\n 			--saveInitBaMMs \n"
			"				Write initialized BaMM(s) to disk.\n\n");
	printf("\n 			--saveBaMMs\n"
			"				Write optimized BaMM(s) to disk.\n\n");
	printf("\n 			-h, --help\n"
			"				 Printout this help function.\n\n");
	printf("\n===============================================================================\n");
}

void Global::destruct(){
    Alphabet::destruct();
    if( alphabetType ) 			delete[] alphabetType;
    if( posSequenceBasename ) 	free( posSequenceBasename );
    if( negSequenceBasename ) 	free( negSequenceBasename );
    if( posSequenceSet )	 	delete posSequenceSet;
    if( negSequenceSet ) 		delete negSequenceSet;
    if( initialModelBasename ) 	free( initialModelBasename );
}

void Global::debug(){

    // check Global Parameter settings:
    fprintf( stdout, "outputDirectory        = %s \n", outputDirectory);
    fprintf( stdout, "posSequenceFilename    = %s \n", posSequenceFilename);
    fprintf( stdout, "posSequenceBasename    = %s \n", posSequenceBasename);
    fprintf( stdout, "negSequenceFilename    = %s \n", negSequenceFilename);
    fprintf( stdout, "negSequenceBasename    = %s \n", negSequenceBasename);
    fprintf( stdout, "\n");
    fprintf( stdout, "negGiven               = %d \n", negGiven);
    fprintf( stdout, "intensityFilename      = %s \n", intensityFilename);
    fprintf( stdout, "alphabetType           = %s \n", alphabetType);
    fprintf( stdout, "revcomp                = %d \n", revcomp);
    fprintf( stdout, "\n");
    fprintf( stdout, "BaMMpatternFilename    = %s \n", BaMMpatternFilename);
    fprintf( stdout, "bindingSiteFilename    = %s \n", bindingSiteFilename);
    fprintf( stdout, "PWMFilename            = %s \n", PWMFilename);
    fprintf( stdout, "BaMMFilename           = %s \n", BaMMFilename);
    fprintf( stdout, "initialModelBasename   = %s \n", initialModelBasename);
    fprintf( stdout, "\n");
    fprintf( stdout, "modelOrder             = %d \n", modelOrder);
    fprintf( stdout, "modelBeta              = %f \n", modelBeta);
    fprintf( stdout, "modelGamma             = %f \n", modelGamma);
    fprintf( stdout, "modelAlpha             =");
    for( int k = 0; k < modelOrder + 1; k++ ){
        fprintf( stdout, " %f", modelAlpha[k]);
    }
    fprintf( stdout, " \n");
    fprintf( stdout, "addColumns             = %d %d \n", addColumns.at(0), addColumns.at(1));
    fprintf( stdout, " \n");
    fprintf( stdout, "bgModelAlpha           =");
    for( int k = 0; k < bgModelOrder + 1; k++ ){
        fprintf( stdout, " %f", bgModelAlpha[k]);
    }
    fprintf( stdout, " \n");
    fprintf( stdout, "bgModelOrder           = %d \n", bgModelOrder);
    fprintf( stdout, "maxEMIterations        = %d \n", maxEMIterations);
    fprintf( stdout, "epsilon                = %f \n", epsilon);
    fprintf( stdout, "noAlphaOptimization    = %d \n", noAlphaOptimization);
    fprintf( stdout, "noQOptimization        = %d \n", noQOptimization);
    fprintf( stdout, "setSlow                = %d \n", setSlow);
    fprintf( stdout, "\n");
    fprintf( stdout, "FDR                    = %d \n", FDR);
    fprintf( stdout, "mFold                  = %d \n", mFold);
    fprintf( stdout, "cvFold                 = %d \n", cvFold);
    fprintf( stdout, "\n");
    fprintf( stdout, "verbose                = %d \n", verbose);
    fprintf( stdout, "debugMode              = %d \n", debugMode);
    fprintf( stdout, "saveInitBaMMs          = %d \n", saveInitBaMMs);
    fprintf( stdout, "saveBaMMs              = %d \n", saveBaMMs);
    fprintf( stdout, "\n");
    fprintf( stdout, "posSequenceSet::sequenceFilepath_  = %s \n", posSequenceSet->getSequenceFilepath().c_str());
    fprintf( stdout, "posSequenceSet::intensityFilepath_ = %s \n", posSequenceSet->getIntensityFilepath().c_str());
    fprintf( stdout, "posSequenceSet::N_                 = %d \n", posSequenceSet->getN());
    fprintf( stdout, "posSequenceSet::minL_              = %d \n", posSequenceSet->getMinL());
    fprintf( stdout, "posSequenceSet::maxL_              = %d \n", posSequenceSet->getMaxL());

    fprintf( stdout, "\n\n");
    fprintf( stdout, "negSequenceSet::sequenceFilepath_  = %s \n", negSequenceSet->getSequenceFilepath().c_str());
    fprintf( stdout, "negSequenceSet::intensityFilepath_ = %s \n", negSequenceSet->getIntensityFilepath().c_str());
    fprintf( stdout, "negSequenceSet::N_                 = %d \n", negSequenceSet->getN());
    fprintf( stdout, "negSequenceSet::minL_              = %d \n", negSequenceSet->getMinL());
    fprintf( stdout, "negSequenceSet::maxL_              = %d \n", negSequenceSet->getMaxL());

    fprintf( stdout, "\n\n");
    fprintf( stdout, "posFoldIndices \n");
    for( size_t fold = 0; fold < posFoldIndices.size(); fold++ ){
        fprintf( stdout, "               %d . fold  = ", static_cast<int>( fold ) );
        for( int n = 0; n < std::min( 10,  static_cast<int>( posFoldIndices[fold].size() )); n++ ){
            fprintf( stdout, "%d ", posFoldIndices[fold][n] );
        }
        fprintf( stdout, " ..... ( L = %d )\n", static_cast<int>( posFoldIndices[fold].size() ));
    }
    if( negGiven ){
        fprintf( stdout, "\n");
        fprintf( stdout, "negFoldIndices \n");
        for( size_t fold = 0; fold < negFoldIndices.size() ; fold++ ){
            fprintf( stdout, "               %d . fold  = ", static_cast<int>( fold ) );
            for( int n = 0; n < std::min( 10,  static_cast<int>( negFoldIndices[fold].size() )); n++ ){
                fprintf( stdout, "%d ", negFoldIndices[fold][n] );
            }
            fprintf( stdout, " ..... ( L = %d )\n", static_cast<int>( negFoldIndices[fold].size() ));
        }
    }
    fprintf( stdout, "\n\n");


}
