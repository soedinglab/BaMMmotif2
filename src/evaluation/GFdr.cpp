//
// Created by wanwan on 03.09.17.
//

#include "GFdr.h"

char*               GFdr::outputDirectory = NULL;		// output directory

char*               GFdr::posSequenceFilename = NULL;	// filename of positive sequence FASTA file
std::string         GFdr::posSequenceBasename;			// basename of positive sequence FASTA file
SequenceSet*        GFdr::posSequenceSet = NULL;		// positive sequence set
float               GFdr::q = 0.9f;						// prior probability for a positive sequence to contain a motif
bool                GFdr::ss = false;					// only search on single strand sequences

char*               GFdr::negSequenceFilename = NULL;	// filename of negative sequence FASTA file
SequenceSet*        GFdr::negSequenceSet = NULL;		// negative sequence set

// alphabet options
char*			    GFdr::alphabetType = NULL;			// alphabet type is defaulted to standard which is ACGT

// initial model(s) options
char*			    GFdr::initialModelFilename = NULL; 	// filename of initial model
std::string         GFdr::initialModelTag;				// tag for initializing the model
size_t              GFdr::num = std::numeric_limits<size_t>::max(); // number of init that are to be optimized
bool                GFdr::mops = false;					// learn MOPS model
bool                GFdr::zoops = true;					// learn ZOOPS model

// model options
size_t              GFdr::modelOrder = 2;				// model order
std::vector<float>  GFdr::modelAlpha( modelOrder+1, 1.f );// initial alphas
float               GFdr::modelBeta = 7.0f;				// alpha_k = beta x gamma^k for k > 0
float               GFdr::modelGamma = 3.0f;
std::vector<size_t> GFdr::addColumns = {0,0};			// add columns to the left and right of initial model

bool                GFdr::interpolateBG = true;			// calculate prior probabilities from lower-order probabilities
                                                        // instead of background frequencies of mononucleotides
// background model options
size_t			    GFdr::bgModelOrder = 2;				// background model order, defaults to 2
std::vector<float>  GFdr::bgModelAlpha( bgModelOrder+1, 1.f );// background model alpha

// refinement options
bool                GFdr::EM = false;					// flag to trigger EM learning
bool			    GFdr::CGS = false;					// flag to trigger Collapsed Gibbs sampling

// FDR options
size_t		        GFdr::mFold = 10;					// number of negative sequences as multiple of positive sequences
size_t		        GFdr::cvFold = 5;					// number of cross-validation (cv) folds
size_t		        GFdr::sOrder = 2;					// k-mer order for sampling negative sequence set

// printout options
bool			    GFdr::savePvalues = false;			// write p-values for each log odds score from sequence set
bool			    GFdr::saveLogOdds = false;			// write the log odds of positive and negative sets to disk


void GFdr::init( int nargs, char* args[] ){

    readArguments( nargs, args );

    Alphabet::init( alphabetType );

    // read in positive, negative and background sequence set
    posSequenceSet = new SequenceSet( posSequenceFilename, ss );
    negSequenceSet = new SequenceSet( negSequenceFilename, ss );
}

int GFdr::readArguments( int nargs, char* args[] ){

    /**
     * read in command line to get options process flags from user
     */

    if( nargs < 3 ) {
        std::cerr << "Error: Arguments are missing!" << std::endl;
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
    for( int i = 3; i < nargs; i++ ){
        if( !strcmp( args[i], "-h" ) or !strcmp( args[i], "--help" ) ){
            printHelp();
            exit( 1 );
        } else if( !strcmp( args[i], "-q" ) ){
            if( ++i >= nargs ){
                printHelp();
                std::cerr << "No expression following -q" << std::endl;
                exit( 2 );
            }
            q = std::stof( args[i] );
        } else if( !strcmp( args[i], "--ss" ) ){
            ss = true;
        } else if( !strcmp( args[i], "--negSeqFile" ) ){
            if( ++i >= nargs ){
                printHelp();
                std::cerr << "No expression following --negSeqFile" << std::endl;
                exit( 2 );
            }
            negSequenceFilename = args[i];
        } else if( !strcmp( args[i], "--alphabet" ) ){
            if( ++i >= nargs ){
                printHelp();
                std::cerr << "No expression following --alphabet" << std::endl;
                exit( 2 );
            }
            alphabetType = args[i];
        } else if( !strcmp( args[i], "--bindingSiteFile" ) ){
            if( ++i >= nargs ){
                printHelp();
                std::cerr << "No expression following --bindingSiteFile" << std::endl;
                exit( 2 );
            }
            initialModelFilename = args[i];
            initialModelTag = "bindingsites";
        } else if( !strcmp( args[i], "--PWMFile" ) ){
            if( ++i >= nargs ){
                printHelp();
                std::cerr << "No expression following --PWMFile" << std::endl;
                exit( 2 );
            }
            initialModelFilename = args[i];
            initialModelTag = "PWM";
        } else if( !strcmp( args[i], "--BaMMFile" ) ){
            if( ++i >= nargs ){
                printHelp();
                std::cerr << "No expression following --BaMMFile" << std::endl;
                exit( 2 );
            }
            initialModelFilename = args[i];
            initialModelTag = "BaMM";
        } else if( !strcmp( args[i], "--maxPWM" ) ){
            if( ++i >= nargs ){
                printHelp();
                std::cerr << "No expression following --maxPWM" << std::endl;
                exit( 2 );
            }
            num = std::stoi( args[i] );
        } else if( !strcmp( args[i], "--mops" ) ){
            mops = true;
        } else if( !strcmp( args[i], "-k" ) or !strcmp( args[i], "--order" ) ){
            if( ++i >= nargs ){
                printHelp();
                std::cerr << "No expression following -k/--order" << std::endl;
                exit( 2 );
            }
            modelOrder = std::stoi( args[i] );
        } else if( !strcmp( args[i], "-b" ) or !strcmp( args[i], "--beta" ) ){
            if( ++i >= nargs ){
                printHelp();
                std::cerr << "No expression following -b/--beta" << std::endl;
                exit( 2 );
            }
            modelBeta = std::stof( args[i] );
        } else if( !strcmp( args[i], "-r" ) or !strcmp( args[i], "--gamma" ) ){
            if( ++i >= nargs ){
                printHelp();
                std::cerr << "No expression following -r/--gamma" << std::endl;
                exit( 2 );
            }
            modelGamma = std::stof( args[i] );
        } else if( !strcmp( args[i], "-K" ) or !strcmp( args[i], "-Order" ) ){
            if( ++i >= nargs ){
                printHelp();
                std::cerr << "No expression following -K/--Order" << std::endl;
                exit( 2 );
            }
            bgModelOrder = std::stoi( args[i] );
        } else if( !strcmp( args[i], "--extend" ) ){
            if( ++i >= nargs ){
                printHelp();
                std::cerr << "No expression following --extend" << std::endl;
                exit( 2 );
            }
            addColumns.at(0) = std::stoi( args[i] );
            addColumns.at(1) = std::stoi( args[i] );
            if( ++i < nargs ){
                addColumns.at(1) = std::stoi( args[i] );
            }
        } else if( !strcmp( args[i], "--EM" ) ){
            EM = true;
        } else if( !strcmp( args[i], "--CGS" ) ){
            CGS = true;
        } else if( !strcmp( args[i], "-m" ) or !strcmp( args[i], "--mFold" ) ){
            if( ++i >= nargs ){
                printHelp();
                std::cerr << "No expression following -m/--mFold" << std::endl;
                exit( 2 );
            }
            mFold = std::stoi( args[i] );
        } else if( !strcmp( args[i], "--cvFold" ) ){
            if( ++i >= nargs ){
                printHelp();
                std::cerr << "No expression following --cvFold" << std::endl;
                exit( 2 );
            }
            cvFold = std::stoi( args[i] );
        } else if( !strcmp( args[i], "-s" ) or !strcmp( args[i], "--sOrder" ) ){
            if( ++i >= nargs ){
                printHelp();
                std::cerr << "No expression following -s/--sOrder" << std::endl;
                exit( 2 );
            }
            sOrder = std::stoi( args[i] );
        } else if( !strcmp( args[i], "--savePvalues" ) ){
            savePvalues = true;
        } else if( !strcmp( args[i], "--saveLogOdds" ) ){
            saveLogOdds = true;
        } else {
            std::cerr << "Ignoring unknown option " << args[i] << std::endl;
        }
    }

    if( negSequenceFilename == NULL ){
        negSequenceFilename = posSequenceFilename;
    }

    if( alphabetType == NULL ){
        alphabetType = new char[9];
        strcpy( alphabetType, "STANDARD" );
    }

    if( initialModelFilename == NULL ){
        std::cerr << "Error: No initial model is provided." << std::endl;
        exit( 1 );
    }

    modelAlpha.resize( modelOrder+1 );
    if( modelOrder > 0 ){
        for( size_t k = 1; k < modelOrder+1; k++ ){
            // alpha = beta * gamma^k
            modelAlpha[k] = modelBeta * powf( modelGamma, ( float )k );
        }
    }

    bgModelAlpha.resize( bgModelOrder+1 );
    if( bgModelOrder > 0 ){
        for( size_t k = 1; k < bgModelOrder+1; k++ ){
            bgModelAlpha[k] = 10.0f;
        }
    }

    return 0;
}

void GFdr::printHelp(){
    printf("\n==================================================================\n");
    printf("\n SYNOPSIS:	FDR OUTDIR SEQFILE [options] \n\n");
    printf("\t DESCRIPTION \n");
    printf("		Evaluate models optimized by BaMM and output precision-recall .\n\n");
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
    printf("\n 			-q <FLOAT> \n"
                   "				Prior probability for a positive sequence to contain\n"
                   "				a motif. The default value is 0.9.\n\n");
    printf("\n			--negSeqFile \n"
                   "				FASTA file with negative/background sequences used\n"
                   "				to learn the (homogeneous) background BaMM.\n"
                   "				If not specified, the background BaMM is learned\n"
                   "				from the positive sequences. \n\n");
    printf("\n		Options for initialize BaMM(s) from file: \n");
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
    printf("\n 		Options for optimization: \n");
    printf("\n 			--EM  \n"
                   "				Triggers Expectation Maximization (EM) algorithm.\n "
                   "				Defaults to false.\n\n");
    printf("\n 			--CGS\n"
                   "				Triggers Collapsed Gibbs Sampling (CGS) algorithm.\n"
                   "				Defaults to false.\n\n");
    printf("\n 		Options for output:	\n");
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
    printf("\n 			-h, --help\n"
                   "				Printout this help function.\n\n");
    printf("\n==================================================================\n");
}

void GFdr::destruct(){
    Alphabet::destruct();
    if( alphabetType ) 			delete[] alphabetType;
    if( posSequenceSet )	 	delete posSequenceSet;
    if( negSequenceSet ) 		delete negSequenceSet;
}
