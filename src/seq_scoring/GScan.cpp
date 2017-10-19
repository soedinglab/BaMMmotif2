//
// Created by wanwan on 04.09.17.
//

#include "GScan.h"

char*               GScan::outputDirectory = NULL;		// output directory

char*               GScan::posSequenceFilename = NULL;	// filename of positive sequence FASTA file
std::string         GScan::posSequenceBasename;			// basename of positive sequence FASTA file
SequenceSet*        GScan::posSequenceSet = NULL;		// positive sequence set
float               GScan::q = 0.9f;					// prior probability for a positive sequence to contain a motif
bool                GScan::ss = false;					// only search on single strand sequences

char*               GScan::negSequenceFilename = NULL;	// filename of negative sequence FASTA file
SequenceSet*        GScan::negSequenceSet = NULL;		// negative sequence set

// alphabet options
char*			    GScan::alphabetType = NULL;			// alphabet type is defaulted to standard which is ACGT

// initial model(s) options
char*			    GScan::initialModelFilename = NULL; // filename of initial model
std::string         GScan::initialModelTag;				// tag for initializing the model
size_t              GScan::maxPWM = std::numeric_limits<size_t>::max(); // number of init that are to be optimized

// model options
size_t              GScan::modelOrder = 2;				// model order
std::vector<float>  GScan::modelAlpha( modelOrder+1, 1.f );// initial alphas
float               GScan::modelBeta = 7.0f;			// alpha_k = beta x gamma^k for k > 0
float               GScan::modelGamma = 3.0f;
std::vector<size_t> GScan::addColumns = {0,0};			// add columns to the left and right of initial model

bool                GScan::interpolateBG = true;		// calculate prior probabilities from lower-order probabilities
                                                        // instead of background frequencies of mononucleotides
// background model options
char*				GScan::bgModelFilename = NULL;		// path to the background model file
size_t			    GScan::bgModelOrder = 2;			// background model order, defaults to 2
std::vector<float>  GScan::bgModelAlpha( bgModelOrder+1, 1.f );// background model alpha

float               GScan::pvalCutoff = 0.0001f;        // cutoff of p-value

void GScan::init( int nargs, char* args[] ){

    readArguments( nargs, args );

    Alphabet::init( alphabetType );

    // read in positive, negative and background sequence set
    posSequenceSet = new SequenceSet( posSequenceFilename, ss );
    negSequenceSet = new SequenceSet( negSequenceFilename, ss );
}

int GScan::readArguments( int nargs, char* args[] ){

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
        } else if( !strcmp( args[i], "--bgModelFile" ) ){
            if( ++i >= nargs ){
                printHelp();
                std::cerr << "No expression following --bgModelFile" << std::endl;
                exit( 2 );
            }
            bgModelFilename = args[i];
        } else if( !strcmp( args[i], "--maxPWM" ) ){
            if( ++i >= nargs ){
                printHelp();
                std::cerr << "No expression following --maxPWM" << std::endl;
                exit( 2 );
            }
            maxPWM = std::stoi( args[i] );
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
        } else if( !strcmp( args[i], "-K" ) or !strcmp( args[i], "--Order" ) ){
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
        }  else if( !strcmp( args[i], "--pvalCutoff" ) ){
            if( ++i >= nargs ){
                printHelp();
                std::cerr << "No expression following --pvalCutoff" << std::endl;
                exit( 2 );
            }
            pvalCutoff = std::stof( args[i] );
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
void GScan::printHelp(){
    std::cout << "DESCRIPTION" << std::endl
              << "Scan given sequence set for the query motifs and output motif occurrences" << std::endl << std::endl
              << "SYNOPSIS:\tBaMMScan OUTDIR SEQFILE [options]" << std::endl << std::endl
              << "OPTIONS:" << std::endl
              << "\tOptions for sequence file:" << std::endl
              << "\t\t--ss" << std::endl
              << "\t\t\tSearch motif only on single-strand sequences. " << std::endl
              << "\t\t\tThis option is not recommended for analyzing " << std::endl
              << "\t\t\tChIP-seq data. " << std::endl
              << "\t\t\tBy default, BaMM searches motifs on both strands." << std::endl
              << "\t\t--negSeqFile" << std::endl
              << "\t\t\tFASTA file with negative/background sequences used to " << std::endl
              << "\t\t\tlearn the (homogeneous) background BaMM. " << std::endl
              << "\t\t\tIf not specified, background model is learned from positive sequences." << std::endl
              << "\tOptions for query motifs:" << std::endl
              << "\t\t--bindingSiteFile <STRING>" << std::endl
              << "\t\t\tFile with binding sites of equal length(one per line)." << std::endl
              << "\t\t--PWMFile <STRING>" << std::endl
              << "\t\t\tFile that contains position weight matrices(PWMs)." << std::endl
              << "\t\t\tSame format as MEME" << std::endl
              << "\t\t--BaMMFile <STRING>" << std::endl
              << "\t\t\tFile that contains a model in bamm file format." << std::endl
              << "\t\t--bgModelFile <STRING>" << std::endl
              << "\t\t\tFile that contains a background model in bamm file format." << std::endl
              << "\t\t--maxPWM <INTEGER>" << std::endl
              << "\t\t\tmaximal number of PWMs that should be optimized." << std::endl;
}

void GScan::destruct(){
    Alphabet::destruct();
    if( alphabetType ) 			delete[] alphabetType;
    if( posSequenceSet )	 	delete posSequenceSet;
    if( negSequenceSet ) 		delete negSequenceSet;
}