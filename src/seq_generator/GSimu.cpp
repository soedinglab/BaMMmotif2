//
// Created by wanwan on 18.09.17.
//

#include "GSimu.h"

char*               GSimu::outputDirectory = NULL;		// output directory

char*               GSimu::sequenceFilename = NULL;	    // filename of input sequence FASTA file
std::string         GSimu::sequenceBasename;			// basename of input sequence FASTA file
SequenceSet*        GSimu::sequenceSet = NULL;		    // input sequence set
float               GSimu::q = 0.9f;					// prior probability for a positive sequence to contain a motif

// alphabet options
char*			    GSimu::alphabetType = NULL;			// alphabet type is defaulted to standard which is ACGT

// initial model(s) options
char*			    GSimu::initialModelFilename = NULL; // filename of initial model
std::string         GSimu::initialModelTag;				// tag for initializing the model

// model options
size_t              GSimu::modelOrder = 2;				// model order
std::vector<float>  GSimu::modelAlpha( modelOrder+1, 1.f );// initial alphas
float               GSimu::modelBeta = 7.0f;			// alpha_k = beta x gamma^k for k > 0
float               GSimu::modelGamma = 3.0f;
std::vector<size_t> GSimu::addColumns = {0,0};			// add columns to the left and right of initial model

// instead of background frequencies of mononucleotides
// background model options
size_t			    GSimu::bgModelOrder = 2;			// background model order, defaults to 2
std::vector<float>  GSimu::bgModelAlpha( bgModelOrder+1, 1.f );// background model alpha

// simulation options
size_t		        GSimu::mFold = 10;					// number of negative sequences as multiple of positive sequences
size_t		        GSimu::sOrder = 2;					// k-mer order for sampling negative sequence set
bool                GSimu::sampleBgset = false;
bool                GSimu::maskSeqset = false;
bool                GSimu::embedSeqset = false;
size_t              GSimu::at = 0;                      // default embedding position as 0, later it will be randomized
                                                        // if not chosen by users

void GSimu::init( int nargs, char* args[] ){

    readArguments( nargs, args );

    Alphabet::init( alphabetType );

    // read in the given sequence set
    // only the given single stand
    sequenceSet = new SequenceSet( sequenceFilename, true );

}

int GSimu::readArguments( int nargs, char* args[] ){

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
    sequenceFilename = args[2];
    sequenceBasename = baseName( sequenceFilename );

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
        } else if( !strcmp( args[i], "-m" ) or !strcmp( args[i], "--mFold" ) ){
            if( ++i >= nargs ){
                printHelp();
                std::cerr << "No expression following -m/--mFold" << std::endl;
                exit( 2 );
            }
            mFold = std::stoi( args[i] );
        } else if( !strcmp( args[i], "-s" ) or !strcmp( args[i], "--sOrder" ) ){
            if( ++i >= nargs ){
                printHelp();
                std::cerr << "No expression following -s/--sOrder" << std::endl;
                exit( 2 );
            }
            sOrder = std::stoi( args[i] );
        } else if( !strcmp( args[i], "--sampleBgset" ) ){
            sampleBgset = true;
        } else if( !strcmp( args[i], "--maskSeqset" ) ){
            maskSeqset = true;
        } else if( !strcmp( args[i], "--embedSeqset" ) ){
            embedSeqset = true;
        } else if( !strcmp( args[i], "--at" ) ){
            if( ++i >= nargs ){
                printHelp();
                std::cerr << "No expression following --at" << std::endl;
                exit( 2 );
            }
            at = std::stoi( args[i] );
        } else {
            std::cerr << "Ignoring unknown option " << args[i] << std::endl;
        }
    }

    if( alphabetType == NULL ){
        alphabetType = new char[9];
        strcpy( alphabetType, "STANDARD" );
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

void GSimu::printHelp(){
    printf("\n==================================================================\n");
    printf("\n SYNOPSIS:	BaMMSimu OUTDIR SEQFILE [options] \n\n");
    printf("\t DESCRIPTION \n");
    printf("		Simulate different sequence set.\n\n");
    printf("\t OUTDIR:  output directory for all results. \n");
    printf("\t SEQFILE: file with positive sequence set in FASTA format.\n\n");
    printf("\n OPTIONS: \n");
    printf("\n		Options for reading in sequence file: \n");
    printf("\n			--alphabet <STRING> \n"
                   "				STANDARD.		For alphabet type ACGT, by default;\n"
                   "				METHYLC. 		For alphabet type ACGTM;\n"
                   "				HYDROXYMETHYLC.	For alphabet type ACGTH;\n"
                   "				EXTENDED.		For alphabet type ACGTMH.\n\n");
    printf("\n 			-q <FLOAT> \n"
                   "				Prior probability for a positive sequence to contain\n"
                   "				a motif. The default value is 0.9.\n\n");
    printf("\n 		Options for simulation: \n");
    printf("\n 			--sampleBgset\n"
                   "				Sample background sequence set based on s-mer frequencies \n"
                   "                from the input sequence set. Defaults to false.\n\n");
    printf("\n 			--maskSeqset\n"
                   "				Mask the given motif from the input sequence set.\n"
                   "				Defaults to false.\n\n");
    printf("\n 			--embedSeqset\n"
                   "				Embed the given motif into the input sequence set.\n"
                   "				Defaults to false.\n\n");
    printf("\n 			--at <INTEGER>\n"
                   "				Fix the position for embedding the given motif.\n"
                   "				Defaults to random positions.\n\n");
    printf("\n		Options for initialize BaMM(s) from file: \n");
    printf("\n 			--bindingSiteFile <STRING> \n"
                   "				File with binding sites of equal length(one per line).\n\n");
    printf("\n 			--PWMFile <STRING> \n"
                   "				File that contains position weight matrices(PWMs).\n");
    printf("\n 			--BaMMFile <STRING> \n"
                   "				File that contains a model in bamm file format.\n\n");
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
    printf("\n 			-h, --help\n"
                   "				Printout this help function.\n\n");
    printf("\n==================================================================\n");
}

void GSimu::destruct(){
    Alphabet::destruct();
    if( alphabetType ) 		delete[] alphabetType;
    if( sequenceSet )	 	delete sequenceSet;
}
