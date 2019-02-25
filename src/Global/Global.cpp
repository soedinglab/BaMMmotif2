#include "Global.h"
#ifdef OPENMP
#include <omp.h>
#endif

char*               Global::outputDirectory = NULL;			// output directory
std::string			Global::outputFileBasename;

char*               Global::posSequenceFilename = NULL;		// filename of positive sequence FASTA file
std::string			Global::posSequenceBasename;			// basename of positive sequence FASTA file
SequenceSet*        Global::posSequenceSet = NULL;			// positive sequence set
bool 		        Global::maskPosSequenceSet = false;		// mask motif patterns from positive sequence set
size_t              Global::maxPosN = 1e5;                  // set maximal number of input sequences used for training model
bool                Global::fixedPosN = false;              // flag for using fixed number of input sequences

char*               Global::negSequenceFilename = NULL;		// filename of negative sequence FASTA file
std::string			Global::negSequenceBasename;			// basename of negative sequence FASTA file
SequenceSet*        Global::negSequenceSet = NULL;			// negative sequence set
bool				Global::negSeqGiven = false;			// a flag for the negative sequence given by users
size_t              Global::negN = 5000;                    // number of negative sequences to be generated
size_t              Global::minNegN = 5000;
bool                Global::fixedNegN = false;              // flag for using fixed number of negative sequences
bool                Global::genericNeg = false;             // flag for generating negative sequences based on generic 2nd-bgModel

// weighting options
char*               Global::intensityFilename = NULL;		// filename of intensity file (i.e. for HT-SELEX data)

char*				Global::alphabetType = NULL;			// alphabet type is defaulted to standard which is ACGT
bool                Global::ss = false;						// only search on single strand sequences
std::vector<size_t> Global:: A2powerK;

// initial model(s) options
char*				Global::initialModelFilename = NULL; 	// filename of initial model
std::string			Global::initialModelBasename;			// basename of initial model
std::string			Global::initialModelTag;				// tag for initializing the model
size_t				Global::maxPWM = 100;                   // number of initial models that are to be optimized
bool				Global::mops = false;					// learn MOPS model
bool				Global::zoops = true;					// learn ZOOPS model

// model options
size_t     			Global::modelOrder = 2;					// model order
size_t              Global::maxOrder = 8;                   // maximal model order
std::vector<float> 	Global::modelAlpha( modelOrder+1, 1.f );// initial alphas
float				Global::modelBeta = 7.0f;				// alpha_k = beta x gamma^k for k > 0
float				Global::modelGamma = 3.0f;
std::vector<size_t>	Global::addColumns( 2 );				// add columns to the left and right of initial model
bool                Global::interpolate = true;             // calculate prior probabilities from lower-order probabilities
                                                            // instead of background frequencies of mononucleotides
bool                Global::interpolateBG = true;			// calculate prior probabilities from lower-order probabilities
                                                            // instead of background frequencies of mononucleotides
// background model options
char*				Global::bgModelFilename = NULL;			// path to the background model file
bool				Global::bgModelGiven = false;			// flag to show if the background model is given or not
size_t				Global::bgModelOrder = 2;				// background model order, defaults to 2
std::vector<float>	Global::bgModelAlpha( bgModelOrder+1, 1.f );// background model alpha

// EM options
bool				Global::EM = false;						// flag to trigger EM learning
float               Global::EMepsilon = 0.001f;             // epsilon for EM convergence
size_t              Global::maxIterations = 1000;           // maximal iterations
float				Global::q = 0.3f;						// prior probability for a positive sequence to contain a motif
bool 				Global::optimizeQ = false;				// optimize hyper-parameter q in EM algorithm
float               Global::f = 0.05f;                      // fraction of sequences to be masked
bool                Global::optimizePos = false;            // optimize positional prior in the EM algorithm

// CGS (Collapsed Gibbs sampling) options
bool				Global::CGS = false;					// flag to trigger Collapsed Gibbs sampling
bool				Global::noInitialZ = false;				// enable initializing z with one E-step
bool				Global::noAlphaOptimization = false;	// disable alpha optimization in CGS
bool				Global::GibbsMHalphas = false;			// enable alpha sampling in CGS using Gibbs Metropolis-Hastings
bool				Global::dissampleAlphas = false;		// enable alpha sampling in CGS using discretely sampling
bool				Global::noZSampling = false;			// disable q sampling in CGS
bool				Global::noQSampling = false;			// disable q sampling in CGS
bool				Global::debugAlphas = false;

// FDR options
bool				Global::FDR = false;					// triggers False-Discovery-Rate (FDR) estimation
size_t				Global::mFold = 10;						// number of negative sequences as multiple of positive sequences
size_t				Global::cvFold = 4;						// size of cross-validation folds
size_t				Global::sOrder = 2;						// the k-mer order for sampling negative sequence set

// motif occurrence options
bool                Global::scoreSeqset = false;            // write logOdds Scores of positive sequence set to disk
float 				Global::pvalCutoff = 0.0001f;			// score cutoff for printing log odds scores as motif hit

// sequence simulator options
bool                Global::maskSeqset = false;
bool                Global::sampleBgset = false;
bool                Global::embedSeqset = false;
size_t              Global::at = 0;


// printout options
bool                Global::verbose = false;
bool                Global::debug = false;                  // debug-mode: prints out everything.
bool				Global::saveBaMMs = false;
bool				Global::savePRs = true;					// write the precision, recall, TP and FP
bool				Global::savePvalues = false;			// write p-values for each log odds score from sequence set
bool				Global::saveLogOdds = false;			// write the log odds of positive and negative sets to disk
bool				Global::saveInitialBaMMs = false;		// write out the initial model to disk
bool				Global::saveBgModel = false;			// write out the background model to disk
bool				Global::generatePseudoSet = false;		// test for alpha learning
std::mt19937		Global::rngx;

// flags for developers
bool			    Global::makeMovie = false;              // print out bamms in each iteration while optimizing
bool 				Global::B2 = false;
bool 				Global::B3 = false;
bool 				Global::B3prime = false;
bool                Global::advanceEM = false;
bool                Global::slowEM = false;

// option for openMP
size_t              Global::threads = 4;                    // number of threads to use
bool                Global::parallelizeOverMotifs = false;  // flag for parallelizing over motifs

Global::Global( int nargs, char* args[] ){

    // seed random number
    srand( 42 );
    rngx.seed( 42 );

    readArguments( nargs, args );

    Alphabet::init( alphabetType );

    // initialize A2powerK vector
    for( size_t i = 0; i < modelOrder + maxOrder; i++ ){
        A2powerK.push_back( ipow( Alphabet::getSize(), i ) );
    }

    // read in positive and negative sequence set
    posSequenceSet = new SequenceSet( posSequenceFilename, ss );
    negSequenceSet = new SequenceSet( negSequenceFilename, ss );

    // check if the input sequences are too few
    if( posSequenceSet->getSequences().size() < cvFold ){
        std::cerr << "Error: Input sequences are too few for training! \n" << std::endl;
        exit( 1 );
    }

    // check if the amount of negative sequences are enough
    if( !negSeqGiven ){
        size_t posN = posSequenceSet->getSequences().size();
        bool rest = minNegN % posN;
        if ( posN * mFold < minNegN ) {
            mFold = minNegN / posN + rest;
        }
    }

    // optional: read in sequence intensities (header and intensity columns?)
    if( intensityFilename != 0 ){
        ;// read in sequence intensity
    }

    if( debug ) {
        verbose = true;
    }

    // print out the input parameters
    printPara();
}

int Global::readArguments( int nargs, char* args[] ){

	/**
	 * read command line to get options
	 * process flags from user
	 */

	if( nargs < 3 ) {
		std::cerr << "Error: Arguments are missing! \n" << std::endl;
		printHelp();
		exit( 1 );
	}

	// read in the output directory and create it
	outputDirectory = args[1];
	createDirectory( outputDirectory );

	// read in the positive sequence file
	posSequenceFilename = args[2];
    posSequenceBasename = baseName( posSequenceFilename );

	// read in options from the third argument on
	GetOpt::GetOpt_pp opt( nargs-2, args+2 );

	if( opt >> GetOpt::OptionPresent( 'h', "help" ) ){
		printHelp();
		exit( 1 );
	}

    // read in the basename of output file,
    // if not given, take the basename of input FASTA file
    if( opt >> GetOpt::OptionPresent( "basename" ) ){
        opt >> GetOpt::Option( "basename", outputFileBasename );
    } else {
        outputFileBasename = posSequenceBasename;
    }
    // mask motif patterns from the positive sequence set
    opt >> GetOpt::OptionPresent( "maskPosSequenceSet", maskPosSequenceSet );
    if( opt >> GetOpt::Option( "maxPosN", maxPosN ) ){
        fixedPosN = true;
    }

	// read in negative sequence file
	if( opt >> GetOpt::OptionPresent( "negSeqFile" ) ){
		negSeqGiven = true;
		opt >> GetOpt::Option( "negSeqFile", negSequenceFilename );
	} else {
	    negSequenceFilename = posSequenceFilename;
	}
	negSequenceBasename = baseName( negSequenceFilename );

    if( opt >> GetOpt::Option( "negN", negN ) ){
        fixedNegN = true;
    }
    opt >> GetOpt::OptionPresent( "genericNeg", genericNeg );

    // Alphabet Type
	if( opt >> GetOpt::OptionPresent( "alphabet" ) ){
		opt >> GetOpt::Option( "alphabet", alphabetType );
	} else {
		alphabetType = new char[9];
		strcpy( alphabetType, "STANDARD" );
	}

	opt >> GetOpt::OptionPresent( "ss", ss );

	// for HT-SELEX data
	opt >> GetOpt::Option( "intensityFile", intensityFilename );

	// get initial model files
	std::string tag;
	if ( opt >> GetOpt::OptionPresent( "bindingSiteFile" ) ){
		opt >> GetOpt::Option( "bindingSiteFile", initialModelFilename );
		initialModelTag = "bindingsites";
	} else if ( opt >> GetOpt::OptionPresent( "PWMFile" ) ){
		opt >> GetOpt::Option( "PWMFile", initialModelFilename );
		initialModelTag = "PWM";
	} else if( opt >> GetOpt::OptionPresent( "BaMMFile" ) ){
		opt >> GetOpt::Option( "BaMMFile", initialModelFilename );
		initialModelTag = "BaMM";
	} else {
		fprintf( stderr, "Error: No initial model is provided.\n" );
		exit( 1 );
	}
	initialModelBasename = baseName( initialModelFilename );

	opt >> GetOpt::Option( "maxPWM", maxPWM );
	opt >> GetOpt::OptionPresent( "mops", mops );

	// model options
	opt >> GetOpt::Option( 'k', "order", modelOrder );

	if( opt >> GetOpt::OptionPresent( 'a', "alpha" ) ){
		modelAlpha.clear();
		opt >> GetOpt::Option( 'a', "alpha", modelAlpha );
		if( modelAlpha.size() != modelOrder+1 ){
			if( modelAlpha.size()  > modelOrder+1 ){
				modelAlpha.resize( modelOrder+1 );
			} else {
				modelAlpha.resize( modelOrder+1, modelAlpha.back() );
			}
		}
	} else {
		if( modelAlpha.size() != modelOrder+1 ){
			if( modelAlpha.size() > modelOrder+1 ){
				modelAlpha.resize( modelOrder+1 );
			} else {
				modelAlpha.resize( modelOrder+1, modelAlpha.back() );
			}
		}
		opt >> GetOpt::Option( 'b', "beta", modelBeta );
		opt >> GetOpt::Option( 'r', "gamma", modelGamma );
		if( modelOrder > 0 ){
			for( size_t k = 1; k < modelOrder+1; k++ ){
				// alpha = beta * gamma^k
				modelAlpha[k] = modelBeta * powf( modelGamma, ( float )k );
			}
		}
	}

	if( opt >> GetOpt::OptionPresent( "extend" ) ){
		addColumns.clear();
		opt >> GetOpt::Option( "extend", addColumns );
		if( addColumns.size() < 1 || addColumns.size() > 2 ){
			fprintf( stderr, "--extend format error.\n" );
			exit( 1 );
		}
		if( addColumns.size() == 1 )
			addColumns.resize( 2, addColumns.back() );
	} else {
		addColumns.at(0) = 0;
		addColumns.at(1) = 0;
	}

	// background model options
	if( opt >> GetOpt::Option( "bgModelFile", bgModelFilename ) ){
		bgModelGiven = true;
	}

	opt >> GetOpt::Option( 'K', "Order", bgModelOrder );

	if( opt >> GetOpt::OptionPresent( 'A', "Alpha" ) ){
		bgModelAlpha.clear();
		opt >> GetOpt::Option( 'A', "Alpha", bgModelAlpha );
		if( bgModelAlpha.size() != bgModelOrder+1 ){
			if( bgModelAlpha.size() > bgModelOrder+1 ){
				bgModelAlpha.resize( bgModelOrder+1 );
			} else {
				bgModelAlpha.resize( bgModelOrder+1, bgModelAlpha.back() );
			}
		}
	} else {
		if( bgModelAlpha.size() != bgModelOrder+1 ){
			if( bgModelAlpha.size() > bgModelOrder+1 ){
				bgModelAlpha.resize( bgModelOrder+1 );
			} else {
				bgModelAlpha.resize( bgModelOrder+1, bgModelAlpha.back() );
			}
		}
		if( bgModelOrder > 0 ){
			for( size_t k = 1; k < bgModelOrder+1; k++ ){
				bgModelAlpha[k] = 10.0f;
			}
		}
	}

	// EM options
	opt >> GetOpt::OptionPresent( "EM", EM );
    opt >> GetOpt::Option( "EMepsilon", EMepsilon );
    opt >> GetOpt::Option( "maxIterations", maxIterations );

	// CGS options
	if( opt >> GetOpt::OptionPresent( "CGS", CGS ) ){
		opt >> GetOpt::OptionPresent( "noInitialZ", noInitialZ );
		opt >> GetOpt::OptionPresent( "noAlphaOpti", noAlphaOptimization );
		opt >> GetOpt::OptionPresent( "GibbsMH", GibbsMHalphas );
		opt >> GetOpt::OptionPresent( "dissample", dissampleAlphas );
		opt >> GetOpt::OptionPresent( "noZSampling", noZSampling );
		opt >> GetOpt::OptionPresent( "noQSampling", noQSampling );
	}
	opt >> GetOpt::OptionPresent( "debugAlphas", debugAlphas );
	opt >> GetOpt::OptionPresent( "generatePseudoSet", generatePseudoSet );

	// saturation options
	opt >> GetOpt::Option( 'q', q );

    // masking options
	opt >> GetOpt::Option( 'f', f );

	// FDR options
	if( opt >> GetOpt::OptionPresent( "FDR", FDR ) ){
		opt >> GetOpt::Option( 'm', "mFold", mFold );
		opt >> GetOpt::Option( 'n', "cvFold", cvFold );
		opt >> GetOpt::Option( 's', "sOrder", sOrder );
	}
	// motif occurrence option
	opt >> GetOpt::OptionPresent( "scoreSeqset", scoreSeqset );
	opt >> GetOpt::Option( "pvalCutoff", pvalCutoff );

    // sequence simulator options
    opt >> GetOpt::OptionPresent( "maskSeqset", maskSeqset );
    opt >> GetOpt::OptionPresent( "sampleBgset", sampleBgset );
    if( opt >> GetOpt::OptionPresent( "embedSeqset", embedSeqset ) ){
        opt >> GetOpt::Option( "at", at );
    }

    // printout options
	opt >> GetOpt::OptionPresent( "verbose", verbose );
	opt >> GetOpt::OptionPresent( "debug", debug );

	opt >> GetOpt::OptionPresent( "saveBaMMs", saveBaMMs );
	opt >> GetOpt::OptionPresent( "saveInitialBaMMs", saveInitialBaMMs );
	opt >> GetOpt::OptionPresent( "savePvalues", savePvalues );
	opt >> GetOpt::OptionPresent( "saveLogOdds", saveLogOdds );
	opt >> GetOpt::OptionPresent( "saveBgModel", saveBgModel );

    // flags for developers
    opt >> GetOpt::OptionPresent( "makeMovie", makeMovie );
	opt >> GetOpt::OptionPresent( "optimizeQ", optimizeQ );
    opt >> GetOpt::OptionPresent( "optimizePos", optimizePos );
	opt >> GetOpt::OptionPresent( "B2", B2 );
	opt >> GetOpt::OptionPresent( "B3", B3 );
	opt >> GetOpt::OptionPresent( "B3prime", B3prime );
    opt >> GetOpt::OptionPresent( "advanceEM", advanceEM );
    opt >> GetOpt::OptionPresent( "slowEM", slowEM );

    // option for openMP
    opt >> GetOpt::Option( "threads", threads );
    opt >> GetOpt::OptionPresent( "parallizeOverMotifs", parallelizeOverMotifs );

#ifdef OPENMP
    omp_set_num_threads( threads );
#endif

	// for remaining unknown options
	if( opt.options_remain() ){
		printHelp();
		std::cerr << "Oops! Unknown option(s) remaining... \n\n";
		exit( 1 );
	}

	return 0;
}

void Global::printPara(){

    std::cout << " ____________________" << std::endl
              << "|*                  *|" << std::endl
              << "| Parameter Settings |" << std::endl
              << "|*__________________*|" << std::endl
              << std::endl;

    std::cout << "Alphabet type\t\t\t\t" << Alphabet::getAlphabet() << std::endl;

    // for positive sequence set
    std::cout << "Positive sequence set\t\t\t"
              << posSequenceFilename << std::endl
              << "Sequence counts\t\t\t\t"
              << posSequenceSet->getSequences().size()
              << " (max.length: "   << posSequenceSet->getMaxL()
              << ", min.length: "   << posSequenceSet->getMinL()
              << ", total bases: "  << posSequenceSet->getBaseSum()
              << ")" << std::endl
              << "Base frequencies\t\t\t";
    for( size_t i = 0; i < Alphabet::getSize(); i++ ){
        std::cout << posSequenceSet->getBaseFrequencies()[i]
                  << "(" << Alphabet::getAlphabet()[i] << ")"
                  << ' ';
    }
    std::cout << std::endl << std::endl;

    if( advanceEM ){
        std::cout << "\n    " << f * 100
                  << "% of the sequences are used for EM after masking."
                  << std::endl;
    }

    // for negative sequence set
    if( negSeqGiven ){
        std::cout << "Negative sequence set\t\t\t"
                  << negSequenceBasename << std::endl
                  << "Sequence counts\t\t\t\t"
                  << negSequenceSet->getSequences().size()
                  << " (max.length: "   << negSequenceSet->getMaxL()
                  << ", min.length: "   << negSequenceSet->getMinL()
                  << ", total bases: "  << negSequenceSet->getBaseSum()
                  << ")" << std::endl
                  << "Base frequencies\t\t\t";
        for( size_t i = 0; i < Alphabet::getSize(); i++ ){
            std::cout << negSequenceSet->getBaseFrequencies()[i]
                      << "(" << Alphabet::getAlphabet()[i] << ")"
                      << ' ';
        }
        std::cout << std::endl << std::endl;
    } else {
        std::cout << "Negative sequence set\t\t\t"
                  << posSequenceSet->getSequences().size() << " x "<< mFold
                  << " sequences are generated based on cond. prob. of "
                  << sOrder+1 << "-mers in positive set"
                  << std::endl << std::endl;
    }

    // for initial model
    std::cout << "Given initial model\t\t\t"    << initialModelBasename << std::endl
              << "Given model type\t\t\t"       << initialModelTag << std::endl << std::endl
              << "MODEL PARAMETERS"             << std::endl
              << "BaMM model order\t\t\t"       << modelOrder << std::endl
              << "Background model order\t\t\t" << bgModelOrder << std::endl
              << "BaMM is trained on\t\t\t";
    if( ss ){
        std::cout << "single-stranded sequences" << std::endl << std::endl;
    } else{
        std::cout << "double-stranded sequences" << std::endl << std::endl;
    }


    // for further functionalities
    std::cout << "Cross-Validation\t\t\t";
    if( FDR ){
        std::cout << "True ("<< cvFold << " folds)" << std::endl;
    } else {
        std::cout << "False" << std::endl;
    }

    std::cout << "Scan positive sequence set\t\t";
    if( scoreSeqset ){
        std::cout << "True" << std::endl;
    } else {
        std::cout << "False" << std::endl;
    }

}

void Global::printStat(){

    std::cout << " ____________" << std::endl
              << "|*          *|" << std::endl
              << "| Statistics |" << std::endl
              << "|*__________*|" << std::endl
              << std::endl;

}

void Global::printHelp(){
	printf("\n==================================================================\n");
	printf("\n SYNOPSIS:	BaMMmotif OUTDIR SEQFILE [options] \n\n");
	printf("\t DESCRIPTION \n");
	printf("		Learn Bayesian inhomogeneous Markov init(BaMMs) from\n"
			"		high-throughput sequencing data.\n\n");
	printf("\t OUTDIR:  output directory for all results. \n");
	printf("\t SEQFILE: file with positive sequence set in FASTA format.\n\n");
	printf("\n OPTIONS: \n");
	printf("\n		Options for reading in sequence file: \n");
	printf("\n			--alphabet <STRING> \n"
			"				STANDARD.		For alphabet type ACGT, by default;\n"
			"				METHYLC. 		For alphabet type ACGTM;\n"
			"				HYDROXYMETHYLC.	For alphabet type ACGTH;\n"
			"				EXTENDED.		For alphabet type ACGTMH.\n\n");
	printf("\n			--ss \n"
			"				Search motif only on single-strand sequences.\n"
			"				This option is not recommended for analyzing\n"
			"				ChIP-seq data. \n"
			"				By default, BaMM searches motifs on both strands.\n\n");
	printf("\n			--negSeqFile \n"
			"				FASTA file with negative/background sequences used\n"
			"				to learn the (homogeneous) background BaMM.\n"
			"				If not specified, the background BaMM is learned\n"
			"				from the positive sequences. \n\n");
	printf("\n		Options for HT-SELEX data: \n");
	printf("\n			--intensityFile	<STRING> \n"
			"				Intensity file name. (Not implemented yet.) \n\n");
	printf("\n		Options for initialize BaMM(s) from file: \n");
	printf("\n 			--BaMMpatternFile <STRING> \n"
			"				File with IUPAC patterns.(Not implemented yet.) \n\n");
	printf("\n 			--bindingSiteFile <STRING> \n"
			"				File with binding sites of equal length(one per line).\n\n");
	printf("\n 			--PWMFile <STRING> \n"
			"				File that contains position weight matrices(PWMs).\n");
	printf("\n 			--BaMMFile <STRING> \n"
			"				File that contains a model in bamm file format.\n\n");
	printf("\n 			--maxPWM <INTEGER> \n"
			"				Number of init to be learned by BaMM!motif, \n"
			"				specific for PWMs. \n"
			"				By default, all the motifs will be optimized.\n\n");
	printf("\n 			--mops \n"
			"				Learn more-than-one-motif-per-sequence (MOPS) model.\n"
			"				By default, it is set as false.\n\n");
	printf("\n 			--zoops \n"
			"				Learn zero-or-one-motif-per-sequence (ZOOPS) model.\n"
			"				By default, it is set as true.\n\n");
	printf("\n 		Options for the (inhomogeneous) motif BaMM: \n");
	printf("\n 			-k, --order <INTEGER> \n"
			"				Model Order. The default is 2. \n\n");
	printf("\n 			-a, --alpha <FLOAT> [<FLOAT>...] \n"
			"				Order-specific prior strength. The default is 1.0 \n"
			"				(for k = 0) and beta x gamma^k (for k > 0). \n"
			"				The options -b and -r are ignored.\n\n");
	printf("\n 			-b, --beta <FLOAT> \n"
			"				For calculating alphas: beta x gamma^k (for k > 0).\n"
			"				The default is 7.0 (for k > 0) \n");
	printf("\n 			-r, --gamma <FLOAT> \n"
			"				For calculating alphas: beta x gamma^k (for k > 0).\n"
			"				The default is 3.0 (for k > 0) \n");
	printf("\n 			--extend <INTEGER>{1, 2} \n"
			"				Extend BaMMs by adding uniformly initialized positions\n"
			"				to the left and/or right of initial BaMMs.\n "
			"				e.g. invoking with --extend 0 2 adds two positions\n"
			"				to the right of initial BaMMs.\n"
			"				Invoking with --extend 2 adds two positions to both\n"
			"				sides of initial BaMMs.\n"
			"				By default, BaMMs are not being extended.\n\n");
	printf("\n 		Options for the (homogeneous) background BaMM: \n");
	printf("\n 			-K, --Order <INTEGER> \n"
			"				Order. The default is 2.\n"
			"				Order of background model should not exceed order of\n"
			"				motif model.\n\n");
	printf("\n 			-A, --Alpha <FLOAT> \n"
			"				Prior strength. The default value is 10.0.\n\n");
	printf("\n 			--bgModelFile <STRING> \n"
			"				Read in background model from a bamm-formatted file.\n"
			"				Defaults to NULL.\n\n");
	printf("\n 		Options for EM: \n");
	printf("\n 			--EM  \n"
			"				Triggers Expectation Maximization (EM) algorithm.\n "
			"				Defaults to false.\n\n");
	printf("\n 			-q <FLOAT> \n"
			"				Prior probability for a positive sequence to contain\n"
			"				a motif. The default value is 0.9.\n\n");
	printf("\n 			-e, --epsilon <FLOAT> \n"
			"				The EM algorithm is deemed to be converged when the\n"
			"				sum over the absolute differences in probabilities\n"
			"				from successive EM rounds is smaller than epsilon.\n"
			"				The default is 0.001.\n\n");
	printf("\n 			--maxIterations <INTEGER> (*) \n"
			"				Limit the number of EM iterations. *For developers.\n\n");
	printf("\n 			--noAlphaOptimization (*) \n"
			"				disable alpha optimization.\n"
			"				Defaults to false. *For developers.\n\n");
	printf("\n 			--noQOptimization (*) \n"
			"				disable q optimization.\n"
			"				Defaults to false. *For developers.\n\n");
	printf("\n 		Options for CGS: \n");
	printf("\n 			--CGS\n"
			"				Triggers Collapsed Gibbs Sampling (CGS) algorithm.\n"
			"				Defaults to false.\n\n");
	printf("\n 			--maxCGSIterations <INTEGER> (*) \n"
			"				Limit the number of CGS iterations. \n"
			"				It should be larger than 5 and defaults to 100.\n\n");
	printf("\n 			--noAlphaSampling (*) \n"
			"				disable alpha sampling.\n"
			"				Defaults to false. *For developers.\n\n");
	printf("\n 			--noQSampling (*) \n"
			"				disable q sampling.\n"
			"				Defaults to false. *For developers.\n\n");
	printf("\n 		Options for FDR: \n");
	printf("\n 			--FDR\n"
			"				Triggers False-Discovery-Rate (FDR) estimation. \n\n");
	printf("\n 			-m, --mFold <INTEGER>\n"
			"				Number of negative sequences as multiple of positive\n"
			"				sequences. The default is 10.\n\n");
	printf("\n 			-n, --cvFold <INTEGER>\n"
			"				Fold number for cross-validation. \n"
			"				The default is 5, which means the training set is\n"
			"				4-fold of the test set.\n\n"
			"			-s, --sOrder <INTERGER>\n"
			"				The order of k-mer for sampling pseudo/negative set.\n"
			"				The default is 2.\n\n");
	printf("\n 		Options for scoring sequence set:\n");
	printf("\n 			--scoreSeqset \n"
			"				Score the sequence set. \n\n");
	printf("\n 			--pvalCutoff \n"
			"				Cutoff of p-value for scoring the sequence set, in order to \n"
			"				find motif occurrences. \n\n");
	printf("\n 		Options for output:	\n");
	printf("\n 			--verbose \n"
			"				Verbose printouts.\n\n");
	printf("\n 			--saveBaMMs\n"
			"				Write optimized BaMM(s) parameters to disk.\n\n");
	printf("\n 			--saveInitialBaMMs \n"
			"				Save the initial BaMM model(s).\n\n");
	printf("\n 			--savePRs\n"
			"				Write true positives(TP), false positives(FP), \n"
			"				FDR and recall values to disk. Defaults to true.\n\n");
	printf("\n 			--savePvalues\n"
			"				Write p-values for plotting area under the \n"
			"				Sensitivity-FDR curve (AUSFC) to disk.\n\n");
	printf("\n 			--saveLogOdds\n"
			"				Write log odds scores from positive and negative \n"
			"				sets to disk.\n"
			"				The default for this log odds score is zero.\n\n");
	printf("\n 			--saveBgModel\n"
			"				Write background model to disk.\n\n");
	printf("\n 			--scoreSeqset\n"
			"				Find the motif occurrences on the sequences due to\n"
			"				log odds scores and write them out.\n\n");
    printf("\n 			--basename\n"
           "				Specify the basename of output files.\n"
           "				By default, all output files have the same basename as \n"
           "				the input positive FASTA file.\n\n");
	printf("\n 			-h, --help\n"
			"				Printout this help function.\n\n");
	printf("\n==================================================================\n");
}

Global::~Global(){
    Alphabet::destruct();
    if( alphabetType ) 			delete[] alphabetType;
    if( posSequenceSet )	 	delete posSequenceSet;
    if( negSequenceSet ) 		delete negSequenceSet;
}
