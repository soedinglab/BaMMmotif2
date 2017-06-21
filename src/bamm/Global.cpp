#include <sys/stat.h>   			// get file status

#include "Global.h"

char*               Global::outputDirectory = NULL;			// output directory

char*               Global::posSequenceFilename = NULL;		// filename of positive sequence FASTA file
std::string			Global::posSequenceBasename;			// basename of positive sequence FASTA file
SequenceSet*        Global::posSequenceSet = NULL;			// positive sequence set
std::vector<std::vector<size_t>> Global::posFoldIndices;		// sequence indices for positive sequence set

char*               Global::negSequenceFilename = NULL;		// filename of negative sequence FASTA file
std::string			Global::negSequenceBasename;			// basename of negative sequence FASTA file
SequenceSet*        Global::negSequenceSet = NULL;			// negative sequence set
std::vector<std::vector<size_t>> Global::negFoldIndices;		// sequence indices for given negative sequence set
bool				Global::negSeqGiven = false;			// a flag for the negative sequence given by users
// weighting options
char*               Global::intensityFilename = NULL;		// filename of intensity file (i.e. for HT-SELEX data)

char*				Global::alphabetType = NULL;			// alphabet type is defaulted to standard which is ACGT
bool                Global::ss = false;						// only search on single strand sequences

// initial model(s) options
char*               Global::BaMMpatternFilename = NULL;		// filename of BaMMpattern file
char*               Global::bindingSiteFilename = NULL;		// filename of binding sites file
char*               Global::PWMFilename = NULL;				// filename of PWM file
char*               Global::BaMMFilename = NULL;			// filename of Markov model (.bamm) file
std::string			Global::initialModelBasename;			// basename of initial model
int					Global::num = 1;						// number of models that are to be optimized
bool				Global::mops = false;					// learn MOPS model
bool				Global::zoops = true;					// learn ZOOPS model

// model options
int        			Global::modelOrder = 2;					// model order
std::vector<float> 	Global::modelAlpha( modelOrder+1, 1.f );// initial alphas
float				Global::modelBeta = 7.0f;				// alpha_k = beta x gamma^k for k > 0
float				Global::modelGamma = 3.0f;
std::vector<int>    Global::addColumns( 2 );				// add columns to the left and right of initial model
bool                Global::interpolate = true;             // calculate prior probabilities from lower-order probabilities
                                                            // instead of background frequencies of mononucleotides
bool                Global::interpolateBG = true;			// calculate prior probabilities from lower-order probabilities
                                                            // instead of background frequencies of mononucleotides
// background model options
char*				Global::bgModelFilename = NULL;			// path to the background model file
bool				Global::bgModelGiven = false;			// flag to show if the background model is given or not
char*				Global::bgSequenceFilename = NULL;		// path to the sequence file where the background model can be learned
bool				Global::bgSeqGiven = false;				// flag to show if the background sequence set is given or not
SequenceSet*        Global::bgSequenceSet = NULL;			// background sequence set
int        			Global::bgModelOrder = 2;				// background model order, defaults to 2
std::vector<float>	Global::bgModelAlpha( bgModelOrder+1, 1.f );// background model alpha

// EM options
bool				Global::EM = false;						// flag to trigger EM learning
int				 	Global::maxEMIterations = std::numeric_limits<int>::max();  // maximum number of iterations
float               Global::epsilon = 0.001f;				// threshold for likelihood convergence parameter
bool                Global::noQOptimization = false;		// disable q optimization
float				Global::q = 0.9f;						// prior probability for a positive sequence to contain a motif

// CGS (Collapsed Gibbs sampling) options
bool				Global::CGS = false;					// flag to trigger Collapsed Gibbs sampling
int 				Global::maxCGSIterations = 100;			// maximum number of iterations for CGS, it should be larger than 5
bool				Global::noInitialZ = false;				// enable initializing z with one E-step
bool				Global::noAlphaOptimization = false;	// disable alpha optimization in CGS
bool				Global::GibbsMHalphas = false;			// enable alpha sampling in CGS using Gibbs Metropolis-Hastings
bool				Global::dissampleAlphas = false;		// enable alpha sampling in CGS using discretely sampling
bool				Global::noZSampling = false;			// disable q sampling in CGS
bool				Global::noQSampling = false;			// disable q sampling in CGS
float				Global::eta = 0.2f;						// learning rate for Gibbs sampling, only for tuning
int					Global::interval = 10;					// interval for sampling z and q, only for tuning
bool				Global::debugAlphas = false;

// FDR options
bool                Global::FDR = false;					// triggers False-Discovery-Rate (FDR) estimation
int        			Global::mFold = 10;						// number of negative sequences as multiple of positive sequences
size_t        		Global::cvFold = 5;						// size of cross-validation folds
int 				Global::sOrder = 2;						// the k-mer order for sampling negative sequence set

// scoring options
float 				Global::scoreCutoff = 0.0;				// score cutoff for printing logodds scores as motif hit

// printout options
bool                Global::verbose = false;
bool                Global::debugMode = false;              // debug-mode: prints out everything.
bool				Global::saveBaMMs = true;
bool				Global::savePRs = true;					// write the precision, recall, TP and FP
bool				Global::savePvalues = false;			// write p-values for each log odds score from sequence set
bool				Global::saveLogOdds = false;			// write the log odds of positive and negative sets to disk
bool				Global::saveInitialBaMMs = false;		// write out the initial model to disk
bool				Global::saveBgModel = false;			// write out the background model to disk
bool                Global::scoreSeqset = false;            // write logOdds Scores of positive sequence set to disk
int					Global::Yk = 10;						// the counts of numbers in Y_ array
bool				Global::generatePseudoSet = false;		// test for alpha learning
std::mt19937		Global::rngx;

void Global::init( int nargs, char* args[] ){

	readArguments( nargs, args );

	Alphabet::init( alphabetType );

	// read in positive, negative and background sequence set
	posSequenceSet = new SequenceSet( posSequenceFilename, ss );
	negSequenceSet = new SequenceSet( negSequenceFilename, ss );
	bgSequenceSet = new SequenceSet( bgSequenceFilename, ss );

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
		strcpy( alphabetType, "STANDARD" );
	}

	opt >> GetOpt::OptionPresent( "ss", ss );

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

	opt >> GetOpt::Option( "num", num );
	opt >> GetOpt::OptionPresent( "mops", mops );
	opt >> GetOpt::Option( "zoops", zoops );

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
				modelAlpha[k] = modelBeta * powf( modelGamma, static_cast<float>( k ) );
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
	if( opt >> GetOpt::Option( "bgModelFile", bgModelFilename ) ){
		bgModelGiven = true;
	}
	// read in background sequence file
	if( opt >> GetOpt::OptionPresent( "bgSeqFile" ) ){
		bgSeqGiven = true;
		opt >> GetOpt::Option( "bgSeqFile", bgSequenceFilename );
	} else {
	    bgSequenceFilename = posSequenceFilename;
	}

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
		opt >> GetOpt::OptionPresent( "noInitialZ", noInitialZ );
		opt >> GetOpt::OptionPresent( "noAlphaOpti", noAlphaOptimization );
		opt >> GetOpt::OptionPresent( "GibbsMH", GibbsMHalphas );
		opt >> GetOpt::OptionPresent( "dissample", dissampleAlphas );
		opt >> GetOpt::OptionPresent( "noZSampling", noZSampling );
		opt >> GetOpt::OptionPresent( "noQSampling", noQSampling );
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
		opt >> GetOpt::Option( 's', "sOrder", sOrder );
	}
	// scoring option
	opt >> GetOpt::Option( "scoreCutoff", scoreCutoff );


	// printout options
	opt >> GetOpt::OptionPresent( "verbose", verbose );
	opt >> GetOpt::OptionPresent( "debug", debugMode );
	opt >> GetOpt::OptionPresent( "saveBaMMs", saveBaMMs );
	opt >> GetOpt::OptionPresent( "saveInitialBaMMs", saveInitialBaMMs );
	opt >> GetOpt::Option( "savePRs", savePRs );
	opt >> GetOpt::OptionPresent( "savePvalues", savePvalues );
	opt >> GetOpt::OptionPresent( "saveLogOdds", saveLogOdds );
	opt >> GetOpt::OptionPresent( "saveBgModel", saveBgModel );
	opt >> GetOpt::OptionPresent( "scoreSeqset", scoreSeqset );

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
			"				STANDARD.		For alphabet type ACGT, default setting; \n"
			"				METHYLC. 		For alphabet type ACGTM; \n"
			"				HYDROXYMETHYLC.	For alphabet type ACGTH; \n"
			"				EXTENDED.		For alphabet type ACGTMH. \n\n");
	printf("\n			--ss \n"
			"				Search motif only on single strand strands (positive sequences). \n"
			"				This option is not recommended for analyzing ChIP-seq data. \n"
			"				By default, BaMM searches motifs on both strands. \n\n");
	printf("\n			--negSeqFile \n"
			"				FASTA file with negative/background sequences used to learn the\n"
			"				(homogeneous) background BaMM. If not specified, the background BaMM\n"
			"				is learned from the positive sequences. \n\n");
	printf("\n		Options for HT-SELEX data: \n");
	printf("\n			--intensityFile	<STRING> \n"
			"				Intensity file name. (Not implemented yet.) \n\n");
	printf("\n		Options for initialize BaMM(s) from file: \n");
	printf("\n 			--BaMMpatternFile <STRING> \n"
			"				File with IUPAC patterns.(Not implemented yet.) \n\n");
	printf("\n 			--bindingSiteFile <STRING> \n"
			"				File with binding sites of equal length (one per line).\n\n");
	printf("\n 			--PWMFile <STRING> \n"
			"				File that contains position weight matrices (PWMs).\n");
	printf("\n 			--BaMMFile <STRING> \n"
			"				File that contains a model in bamm file format.\n\n");
	printf("\n 			--num <INTEGER> \n"
			"				Number of models to be learned by BaMM!motif, specific for PWMs. \n"
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
			"				Order-specific prior strength. The default is 1.0 (for k = 0) and\n"
			"				beta x gamma^k (for k > 0). The options -b and -r are ignored.\n\n");
	printf("\n 			-b, --beta <FLOAT> \n"
			"				For calculating alphas: beta x gamma^k (for k > 0). \n"
			"				The default is 7.0 (for k > 0) \n");
	printf("\n 			-r, --gamma <FLOAT> \n"
			"				For calculating alphas: beta x gamma^k (for k > 0). \n"
			"				The default is 3.0 (for k > 0) \n");
	printf("\n 			--extend <INTEGER>{1, 2} \n"
			"				Extend BaMMs by adding uniformly initialized positions to the left\n"
			"				and/or right of initial BaMMs. e.g. invoking with --extend 0 2 adds\n"
			"				two positions to the right of initial BaMMs. Invoking with --extend 2\n"
			"				adds two positions to both sides of initial BaMMs. By default, BaMMs\n"
			"				are not being extended.\n\n");
	printf("\n 		Options for the (homogeneous) background BaMM: \n");
	printf("\n 			-K, --Order <INTEGER> \n"
			"				Order. The default is 2.\n"
			"				Order of background model should not exceed order of motif model.\n\n");
	printf("\n 			-A, --Alpha <FLOAT> \n"
			"				Prior strength. The default value is 10.0.\n\n");
	printf("\n 			--bgModelFile <STRING> \n"
			"				Read in background model from a bamm-formatted file. Defaults to NULL.\n\n");
	printf("\n 		Options for EM: \n");
	printf("\n 			--EM  \n"
			"				Triggers Expectation Maximization (EM) algorithm. Defaults to false.\n\n");
	printf("\n 			-q <FLOAT> \n"
			"				Prior probability for a positive sequence to contain a motif.\n"
			"				The default value is 0.9.\n\n");
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
			"				Triggers Collapsed Gibbs Sampling (CGS) algorithm. Defaults to false.\n\n");
	printf("\n 			--maxCGSIterations <INTEGER> (*) \n"
			"				Limit the number of CGS iterations. \n"
			"				It should be larger than 5 and defaults to 100.\n\n");
	printf("\n 			--noAlphaSampling (*) \n"
			"				disable alpha sampling. Defaults to false. *For developers.\n\n");
	printf("\n 			--noQSampling (*) \n"
			"				disable q sampling. Defaults to false. *For developers.\n\n");
	printf("\n 		Options for FDR: \n");
	printf("\n 			--FDR\n"
			"				Triggers False-Discovery-Rate (FDR) estimation. \n\n");
	printf("\n 			-m, --mFold <INTEGER>\n"
			"				Number of negative sequences as multiple of positive sequences.\n"
			"				The default is 10.\n\n");
	printf("\n 			-n, --cvFold <INTEGER>\n"
			"				Fold number for cross-validation. \n"
			"				The default is 5, which means the training set is 4-fold of the test set.\n\n"
			"			-s, --sOrder <INTERGER>\n"
			"				The order of k-mer for sampling pseudo/negative set. The default is 2.\n\n");
	printf("\n 		Options for scoring sequence set:\n");
	printf("\n 			--scoreSeqset \n"
			"				Score positive sequence set. \n\n");
	printf("\n 		Options for output:	\n");
	printf("\n 			--verbose \n"
			"				Verbose printouts.\n\n");
	printf("\n 			--saveBaMMs\n"
			"				Write optimized BaMM(s) parameters to disk.\n\n");
	printf("\n 			--saveInitialBaMMs \n"
			"				Save the initial BaMM model(s).\n\n");
	printf("\n 			--savePRs\n"
			"				Write true positives(TP), false positives(FP), FDR and recall values to disk.\n"
			"				Defaults to true.\n\n");
	printf("\n 			--savePvalues\n"
			"				Write p-values for plotting area under the Sensitivity-FDR curve (AUSFC)\n"
			"				to disk.\n\n");
	printf("\n 			--saveLogOdds\n"
			"				Write log odds scores from positive and negative sets to disk.\n\n");
	printf("\n 			--saveBgModel\n"
			"				Write background model to disk.\n\n");
	printf("\n 			-h, --help\n"
			"				Printout this help function.\n\n");
	printf("\n============================================================================================\n");
}

void Global::destruct(){
    Alphabet::destruct();
    if( alphabetType ) 			delete[] alphabetType;
    if( posSequenceSet )	 	delete posSequenceSet;
    if( negSequenceSet ) 		delete negSequenceSet;
    if( bgSequenceSet ) 		delete bgSequenceSet;
}

void Global::debug(){

}
