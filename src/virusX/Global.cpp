#include "Global.h"

char*               Global::inputDirectory = NULL;			// output directory
char*               Global::outputDirectory = NULL;			// output directory

char*               Global::negSequenceFilename = NULL;		// filename of negative sequence FASTA file
char*				Global::negSequenceBasename = NULL;		// basename of negative sequence FASTA file

char*				Global::extension = NULL;				// extension of files in FASTA format, default is fasta

char*				Global::alphabetString = NULL;			// defaults to STANDARD which corresponds to ACGT
bool                Global::revcomp = false;				// also search on reverse complement of sequences

SequenceSet*        Global::negSequenceSet = NULL;			// negative Sequence Set

// background model options
int        			Global::bgModelOrder = 2;				// background model order, defaults to 2
float				Global::bgModelAlpha = 10.0f;			// background model alpha

std::vector< std::vector<int> > Global::negFoldIndices( 1 ); 	// sequence indices for each cross-validation fold

bool                Global::verbose = false;            	// verbose printouts

int* 				Global::powA = NULL;

void Global::init( int nargs, char* args[] ){

	readArguments( nargs, args );

	Alphabet::init( alphabetString );

	// read in or generate negative sequence set
//	if( negSequenceFilename == NULL ){
//		std::cerr << "Error: No sequences provided..." << std::endl;
//		exit( -1 );
//	} else{
//		negSequenceSet = new SequenceSet( negSequenceFilename );
//	}

	// generate folds (fill negFoldIndices)
//	generateFolds();

//	powA = new int[modelOrder+2];
//	for( int k = 0; k < modelOrder+2; k++ )
//		powA[k] = ipow( Alphabet::getSize(), k );
}

int Global::readArguments( int nargs, char* args[] ){

	if( nargs < 2 ) {
		std::cerr << "Error: Less than two arguments provided..." << std::endl;
		printHelp();
		exit( -1 );
	}

	// read in input directory
	inputDirectory = args[1];

	GetOpt::GetOpt_pp opt( nargs, args );

	if( opt >> GetOpt::OptionPresent( 'h', "help" ) ){
		printHelp();
		exit( -1 );
	}

	alphabetString = new char[9];
	strcpy( alphabetString, "STANDARD" );

	if( opt >> GetOpt::OptionPresent( "extension" ) ){
		opt >> GetOpt::Option( "extension", extension );
	} else {
		extension = new char[6];
		strcpy( extension, "fasta");
	}

	DIR *dir;
	struct dirent *ent;
	if( ( dir = opendir ( inputDirectory ) ) != NULL ){
		while( ( ent = readdir ( dir ) ) != NULL ){
			char* ext = strrchr( ent->d_name, '.' );
			if( strcmp( ext+1, extension ) == 0 ){
				printf( "%s\n", ent->d_name );
			}
		}
		closedir( dir );
	} else {
		perror( "" );
		return EXIT_FAILURE;
	}

	// read in the positive sequence file
//	negSequenceFilename = args[2];
//	negSequenceBasename = baseName( negSequenceFilename );

//
//	GetOpt::GetOpt_pp opt( nargs, args );
//
//	// negative sequence set
//	if( opt >> GetOpt::OptionPresent( "negSeqFile" ) ){
//		opt >> GetOpt::Option( "negSeqFile", negSequenceFilename );
//		negSequenceBasename = baseName( negSequenceFilename );
//	}
//
//	// Alphabet Type
//	if( opt >> GetOpt::OptionPresent( "alphabet" ) ){
//		opt >> GetOpt::Option( "alphabet", alphabetType );
//	} else {
//		alphabetType = new char[9];
//		strcpy( alphabetType, "STANDARD");
//	}
//
//	opt >> GetOpt::OptionPresent( "revcomp", revcomp );
//
//	// for HT-SLEX data
//	opt >> GetOpt::Option( "intensityFile", intensityFilename );
//
//	// get initial motif files
//	if( opt >> GetOpt::OptionPresent( "BaMMpatternFile" ) ){
//		opt >> GetOpt::Option( "BaMMpatternFile", BaMMpatternFilename );
//		initialModelBasename = baseName( BaMMpatternFilename );
//	} else if ( opt >> GetOpt::OptionPresent( "bindingSitesFile" ) ){
//		opt >> GetOpt::Option( "bindingSitesFile", bindingSitesFilename );
//		initialModelBasename = baseName( bindingSitesFilename );
//	} else if ( opt >> GetOpt::OptionPresent( "PWMFile" ) ){
//		opt >> GetOpt::Option( "PWMFile", PWMFilename );
//		initialModelBasename = baseName( PWMFilename );
//	} else if( opt >> GetOpt::OptionPresent( "BMMFile" ) ){
//		opt >> GetOpt::Option( "BaMMFile", BaMMFilename );
//		initialModelBasename = baseName( BaMMFilename );
//	} else {
//		fprintf( stderr, "Error: No initial model is provided.\n" );
//		exit( -1 );
//	}
//
//	// model options
//	opt >> GetOpt::Option( 'k', "modelOrder", modelOrder );
//	if( opt >> GetOpt::OptionPresent( "modelAlpha" ) || opt >> GetOpt::OptionPresent( 'a' ) )
//		opt >> GetOpt::Option( 'a', "modelAlpha", modelAlpha );
//	else 		// defaults modelAlpha.at(k) to 10.0
//		for( int k = 0; k < modelOrder+1; k++ )
//			modelAlpha.push_back( 10.0f );
//
//	if( opt >> GetOpt::OptionPresent( "addColumns" ) ){
//		addColumns.clear();
//		opt >> GetOpt::Option( "addColumns", addColumns );
//		if( addColumns.size() < 1 || addColumns.size() > 2 ){
//			fprintf( stderr, "--addColumns format error.\n" );
//			exit( -1 );
//		}
//		if( addColumns.size() == 1 )
//			addColumns.resize( 2, addColumns.back() );
//	} else {
//		addColumns.at(0) = 0;
//		addColumns.at(1) = 0;
//	}
//
//	// background model options
//	opt >> GetOpt::Option( 'K', "bgModelOrder", bgModelOrder );
//	opt >> GetOpt::Option( 'A', "bgModelAlpha", bgModelAlpha );
//
//	// em options
//	opt >> GetOpt::Option( "maxEMIterations", maxEMIterations );
//	opt >> GetOpt::Option( "epsilon", epsilon );
//	opt >> GetOpt::OptionPresent( "noAlphaOptimization", noAlphaOptimization );
//	opt >> GetOpt::OptionPresent( "noQOptimization", noQOptimization );
//
//	// FDR options
//	if( opt >> GetOpt::OptionPresent( "FDR", FDR) ){
//		opt >> GetOpt::Option( "m", mFold  );
//		opt >> GetOpt::Option( "n", cvFold );
//	}
//
//	opt >> GetOpt::OptionPresent( "verbose", verbose );
//
//	opt >> GetOpt::OptionPresent( "setSlow", setSlow );			// to be deleted when release
//
//	if( opt.options_remain() ){
//		std::cerr << "Warning: Unknown options remaining..." << std::endl;
//		return false;
//	}

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

void Global::generateFolds(){
	std::vector<int> v( negSequenceSet->getN() );
	std::iota (v.begin(), v.end(), 0 );

	for( int i = 0; i < negSequenceSet->getN(); i++ ){
		negFoldIndices[0].push_back( i );
	}
}

// calculate the power for integer base
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

	bool developerHelp = true; // print developer-specific options

	printf( "\n" );
	printf( "SYNOPSIS\n" );
	printf( "      VirusX DIRPATH [OPTIONS]\n\n" );
	printf( "DESCRIPTION\n" );
	printf( "      Learn Bayesian homogeneous Markov models from bacterial genomes.\n\n" );
	printf( "      DIRPATH\n"
			"          Input directory with bacterial genomes in FASTA format (one file per\n"
			"	       genome). The default filename extension searched for is <fasta>.\n"
			"          Use option -e (--extension) to change the default setting.\n\n" );
	printf("OPTIONS\n");
	printf( "  Options for input files\n" );
	printf( "      -e|--extension <STRING>\n"
			"          The filename extension of FASTA files searched for in the input\n"
			"          directory. The default is <fasta>.\n\n" );
	printf( "  Options for Bayesian homogeneous Markov models\n" );
	printf( "      -k <INTEGER>\n"
			"          Order. The default is 2.\n\n" );
	printf( "      -a|--alpha <FLOAT>\n"
			"          Prior strength. The default is 10.0.\n\n" );

	if( developerHelp ){
		printf( "  Extend options for Bayesian homogeneous Markov models\n" );
		printf( "      -a|--alpha <FLOAT> [<FLOAT>...] (*)\n"
				"          Order-specific prior strength. The default is 1.0 (for k = 0) and\n"
				"          20 x 3^(k-1) (for k > 0). The options -b and -g are ignored.\n\n" );
		printf( "      -b|--beta <FLOAT> (*)\n"
				"          Calculate order-specific alphas according to beta x gamma^(k-1) (for\n"
				"          k > 0). The default is 20.0.\n\n" );
		printf( "      -g|--gamma <FLOAT> (*)\n"
				"          Calculate order-specific alphas according to beta x gamma^(k-1) (for\n"
				"          k > 0). The default is 3.0.\n\n" );
	}
	printf( "  Output options\n" );
	printf( "      -o|--outputDirectory <DIRPATH>\n"
			"          Output directory for Bayesian homogenous Markov models. The input\n"
			"          directory is the default output directory.\n\n" );
	printf( "      --verbose\n"
			"          Verbose terminal printouts.\n\n" );
	printf( "      -h, --help\n"
			"          Printout this help.\n\n" );
	if( developerHelp ){
		printf( "      (*) Developer options\n\n" );
	}
}

void Global::destruct(){

	Alphabet::destruct();

	if( powA ){
		delete[] powA;
	}
	if( alphabetString ){
		delete[] alphabetString;
	}
	if( negSequenceBasename ){
		free( negSequenceBasename );
	}
}
