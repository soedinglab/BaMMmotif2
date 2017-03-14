#include <sys/stat.h>   			// get file status

#include "Global.h"

char*               Global::outputDirectory = NULL;			// output directory

char*               Global::posSequenceFilename = NULL;		// filename of positive sequence FASTA file
std::string			Global::posSequenceBasename;			// basename of positive sequence FASTA file
SequenceSet*        Global::posSequenceSet = NULL;			// positive Sequence Set
std::vector<std::vector<int>> Global::posFoldIndices;		// sequence indices for positive sequence set

char*               Global::negSequenceFilename = NULL;		// filename of negative sequence FASTA file
std::string			Global::negSequenceBasename;			// basename of negative sequence FASTA file
SequenceSet*        Global::negSequenceSet = NULL;			// negative Sequence Set
std::vector<std::vector<int>> Global::negFoldIndices;		// sequence indices for given negative sequence set
bool				Global::negSeqGiven = false;			// a flag for the negative sequence given by users
// weighting options
char*               Global::intensityFilename = NULL;		// filename of intensity file (i.e. for HT-SELEX data)

char*				Global::alphabetType = NULL;			// alphabet type is defaulted to standard which is ACGT
bool                Global::revcomp = false;				// also search on reverse complement of sequences

// initial model(s) options
char*               Global::BaMMpatternFilename = NULL;		// filename of BaMMpattern file
char*               Global::bindingSiteFilename = NULL;		// filename of binding sites file
char*               Global::PWMFilename = NULL;				// filename of PWM file
char*               Global::BaMMFilename = NULL;			// filename of Markov model (.bamm) file
std::string			Global::initialModelBasename;			// basename of initial model

// model options
int        			Global::modelOrder = 2;					// model order
std::vector<float> 	Global::modelAlpha( modelOrder+1, 1.0f );	// initial alphas
float				Global::modelBeta = 20.0f;				// alpha_k = beta x gamma^k for k > 0
float				Global::modelGamma = 3.0f;
std::vector<int>    Global::addColumns( 2 );				// add columns to the left and right of initial model
bool                Global::interpolate = true;             // calculate prior probabilities from lower-order probabilities
                                                            // instead of background frequencies of mononucleotides
bool                Global::interpolateBG = true;           // calculate prior probabilities from lower-order probabilities
                                                            // instead of background frequencies of mononucleotides
// background model options
int        			Global::bgModelOrder = 2;				// background model order, defaults to 2
std::vector<float>	Global::bgModelAlpha( bgModelOrder+1, 1.0f );	// background model alpha

// EM options
bool				Global::EM = false;						// flag to trigger EM learning
int				 	Global::maxEMIterations = std::numeric_limits<int>::max();  // maximum number of iterations
float               Global::epsilon = 0.001f;				// threshold for likelihood convergence parameter
bool                Global::noQOptimization = false;		// disable q optimization
float				Global::q = 0.9f;						// prior probability for a positive sequence to contain a motif

// CGS (Collapsed Gibbs sampling) options
bool				Global::CGS = false;					// flag to trigger Collapsed Gibbs sampling
int 				Global::maxCGSIterations = 200;			// maximum number of iterations for CGS
bool				Global::noAlphaOptimization = false;	// disable alpha optimization in CGS
bool				Global::alphaSampling = false;			// enable alpha sampling in CGS
bool				Global::noZQSampling = false;			// disable q sampling in CGS
float				Global::eta = 0.2f;						// learning rate for Gibbs sampling, only for tuning
int					Global::interval = 10;					// interval for sampling z and q, only for tuning
bool				Global::debugAlphas = false;


// FDR options
bool                Global::FDR = false;					// triggers False-Discovery-Rate (FDR) estimation
int        			Global::mFold = 10;						// number of negative sequences as multiple of positive sequences
int        			Global::cvFold = 5;						// size of cross-validation folds
int 				Global::sOrder = 2;						// the k-mer order for sampling negative sequence set

// printout options
bool                Global::verbose = false;
bool                Global::debugMode = false;              // debug-mode: prints out everything.
bool				Global::saveBaMMs = true;
bool				Global::savePRs = false;				// write the precision, recall, TP and FP
bool				Global::savePvalues = true;				// write p-values for each log odds score from sequence set
bool				Global::saveLogOdds = false;			// write the log odds of positive and negative sets to disk
bool				Global::saveInitialModel = false;		// write out the initial model to disk
bool				Global::saveBgModel = false;			// write out the background model to disk
int					Global::Yk = 10;						// the counts of numbers in Y_ array
bool				Global::generatePseudoSet = false;		// test for alpha learning

void Global::init( int nargs, char* args[] ){

	readArguments( nargs, args );

	Alphabet::init( alphabetType );

	// read in positive and negative sequence set
	posSequenceSet = new SequenceSet( posSequenceFilename, revcomp );
	negSequenceSet = new SequenceSet( negSequenceFilename, true );

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

	if( nargs < 3 ) {
		std::cerr << "Error: Arguments are missing! \n" << std::endl;
		printHelp();
		exit( -1 );
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
		exit( -1 );
	}

	// read in negative sequence file
	if( opt >> GetOpt::OptionPresent( "negSeqFile" ) ){
		negSeqGiven = true;
		opt >> GetOpt::Option( "negSeqFile", negSequenceFilename );
	} else {
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
	} else if( opt >> GetOpt::OptionPresent( "BaMMFile" ) ){
		opt >> GetOpt::Option( "BaMMFile", BaMMFilename );
		initialModelBasename = baseName( BaMMFilename );
	} else {
		fprintf( stderr, "Error: No initial model is provided.\n" );
		exit( -1 );
	}

	opt >> GetOpt::OptionPresent( "saveInitialModel", saveInitialModel );


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
		opt >> GetOpt::Option( 'b', "beta", modelBeta );
		opt >> GetOpt::Option( 'r', "gamma", modelGamma );
		if( modelOrder > 0 ){
			for( int k = 1; k < modelOrder + 1; k++ ){
				// alpha = beta * gamma^k
				modelAlpha[k] = modelBeta * powf( modelGamma, static_cast<float>( k-1 ) );
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
			}
		}
	}

	// EM options
	if( opt >> GetOpt::OptionPresent( "EM", EM ) ){
		opt >> GetOpt::Option( "maxEMIterations", maxEMIterations );
		opt >> GetOpt::Option( 'e', "epsilon", epsilon );
		opt >> GetOpt::OptionPresent( "noQOptimization", noQOptimization );
	}

	// CGS options
	if( opt >> GetOpt::OptionPresent( "CGS", CGS ) ){
		opt >> GetOpt::Option( "maxCGSIterations", maxCGSIterations );
		opt >> GetOpt::OptionPresent( "noAlphaOptimization", noAlphaOptimization );
		opt >> GetOpt::OptionPresent( "AlphaSampling", alphaSampling );
		opt >> GetOpt::OptionPresent( "noZQSampling", noZQSampling );
		opt >> GetOpt::Option( "eta", eta );
		opt >> GetOpt::Option( "interval", interval );
	}
	opt >> GetOpt::OptionPresent( "debugAlphas", debugAlphas );
	opt >> GetOpt::OptionPresent( "generatePseudoSet", generatePseudoSet );

	// saturation options
	opt >> GetOpt::Option( 'q', q );

	// FDR options
	if( opt >> GetOpt::OptionPresent( "FDR", FDR ) ){
		opt >> GetOpt::Option( 'm', "mFold", mFold  );
		opt >> GetOpt::Option( 'n', "cvFold", cvFold );
		opt >> GetOpt::Option( 's', "samplingOrder", sOrder );
	}

	// printout options
	opt >> GetOpt::OptionPresent( "verbose", verbose );
	opt >> GetOpt::OptionPresent( "debug", debugMode );
	opt >> GetOpt::OptionPresent( "saveBaMMs", saveBaMMs );
	opt >> GetOpt::OptionPresent( "savePRs", savePRs );
	opt >> GetOpt::OptionPresent( "savePvalues", savePvalues );
	opt >> GetOpt::OptionPresent( "saveLogOdds", saveLogOdds );
	opt >> GetOpt::OptionPresent( "saveBgModel", saveBgModel );

	// for remaining unknown options
	if( opt.options_remain() ){
		printHelp();
		std::cerr << "Oops! Unknown option(s) remaining... \n\n";
		exit( -1 );
	}

	return 0;
}

void Global::printHelp(){
	printf("\n============================================================================================\n");
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
	printf("\n			--negSeqFile \n"
			"				negative sequence file. \n\n");
	printf("\n		Options for HT-SELEX data: \n");
	printf("\n			--intensityFile	<STRING> \n"
			"				intensity file name. \n\n");
	printf("\n		Options for initial model: \n");
	printf("\n 			--BaMMpatternFile <STRING> \n"
			"				file that contains patterns.\n\n");
	printf("\n 			--bindingSiteFile <STRING> \n"
			"				file that contains binding sites.\n\n");
	printf("\n 			--PWMFile <STRING> \n"
			"				file that contains PWM data.\n");
	printf("\n 			--BaMMFile <STRING> \n"
			"				file that contains BaMM data.\n\n");
	printf("\n 			--saveInitialModel \n"
			"				save initial model.\n\n");
	printf("\n 		Options for inhomogeneous BaMM: \n");
	printf("\n 			-k, --order <INTEGER> \n"
			"				model Order. The default is 2. \n\n");
	printf("\n 			-a, --alpha <FLOAT> [<FLOAT>...] \n"
			"				Order-specific prior strength. The default is 1.0 (for k = 0) and\n"
			"				beta x gamma^k (for k > 0). The options -b and -r are ignored.\n\n");
	printf("\n 			-b, --beta <FLOAT> \n"
			"				parameter for calculating alphas: beta x gamma^k (for k > 0). \n"
			"				The default is 20 (for k > 0) \n");
	printf("\n 			-r, --gamma <FLOAT> \n"
			"				parameter for calculating alphas: beta x gamma^k (for k > 0). \n"
			"				The default is 3 (for k > 0) \n");
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
	printf("\n 			--EM  \n"
			"				trigger Expectation Maximization (EM) algorithm.\n\n");
	printf("\n 			-q <FLOAT> \n"
			"				Prior probability for a positive sequence to contain a motif.\n"
			"				The default is 0.9.\n\n");
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
	printf("\n 		Options for CGS: \n");
	printf("\n 			--CGS\n"
			"				trigger Collapsed Gibbs Sampling (CGS) algorithm.\n\n");
	printf("\n 			--noAlphaSampling (*) \n"
			"				disable alpha sampling. Defaults to false. *For developers.\n\n");
	printf("\n 			--noQSampling (*) \n"
			"				disable q sampling. Defaults to false. *For developers.\n\n");
	printf("\n 		Options for FDR: \n");
	printf("\n 			--FDR\n"
			"				triggers False-Discovery-Rate (FDR) estimation. The default is false.\n\n");
	printf("\n 			-m, --mFold <INTEGER>\n"
			"				number of negative sequences as multiple of positive sequences.\n"
			"				The default is 10.\n\n");
	printf("\n 			-n, --cvFold <INTEGER>\n"
			"				number of cross-validation folds. The default is 5.\n\n"
			"			-s, --samplingOrder <INTERGER>\n"
			"				order of k-mer for sampling negative set. The default is 2.\n\n");
	printf("\n 		Options for output:	\n");
	printf("\n 			--verbose \n"
			"				verbose printouts. Defaults to false.\n\n");
	printf("\n 			--saveBaMMs\n"
			"				Write optimized BaMM(s) parameters to disk.\n\n");
	printf("\n 			--savePRs\n"
			"				Write true positives(TP), false positives(FP), FDR and recall values to disk.\n\n");
	printf("\n 			--savePvalues\n"
			"				Write p-values for plotting area under the Sensitivity-FDR curve (AUSFC)\n"
			"				to disk.\n\n");
	printf("\n 			--saveLogOdds\n"
			"				Write log odds scores from positive and negative sets to disk.\n\n");
	printf("\n 			--saveBgModel\n"
				"			Write background model to disk.\n\n");
	printf("\n 			-h, --help\n"
			"				Printout this help function.\n\n");
	printf("\n============================================================================================\n");
}

void Global::destruct(){
    Alphabet::destruct();
    if( alphabetType ) 			delete[] alphabetType;
    if( posSequenceSet )	 	delete posSequenceSet;
    if( negSequenceSet ) 		delete negSequenceSet;
}

void Global::debug(){

}
