#include "GPrediction.h"

char*		Global::inputDirectoryBaMMs = NULL;	// input directory with BaMM files
char*		Global::inputDirectorySeqs = NULL;	// input directory with FASTA files

char*		Global::outputDirectory = NULL;		// output directory

char*		Global::extensionBaMMs = NULL;		// extension of BaMM files, defaults to hbcp
char*		Global::extensionSeqs = NULL;		// extension of FASTA files, defaults to fasta

bool		Global::verbose = false;			// verbose printouts

std::string	Global::name="predict";				// program name
std::string	Global::version="0.1.1";			// program version number

void Global::init( int nargs, char* args[] ){

	readArguments( nargs, args );

	char* alphabetString = strdup( "STANDARD" );
	Alphabet::init( alphabetString );
	free( alphabetString );
}

void Global::destruct(){

	if( inputDirectoryBaMMs != NULL ){
		free( inputDirectoryBaMMs );
	}
	if( inputDirectorySeqs != NULL ){
		free( inputDirectorySeqs );
	}
	if( outputDirectory != NULL ){
		free( outputDirectory );
	}
	if( extensionBaMMs != NULL ){
		free( extensionBaMMs );
	}
	if( extensionSeqs != NULL ){
		free( extensionSeqs );
	}

	Alphabet::destruct();
}

int Global::readArguments( int nargs, char* args[] ){

	GetOpt::GetOpt_pp opt( nargs, args );

	if( opt >> GetOpt::OptionPresent( 'h', "help" ) ){
		printHelp();
		exit( -1 );
	}

	if( opt >> GetOpt::OptionPresent( "version" ) ){
		std::cout << name << " version " << version << std::endl;
		exit( -1 );
	}

	if( opt >> GetOpt::OptionPresent( 'e', "extensionBaMMs" ) ){
		opt >> GetOpt::Option( 'e', "extensionBaMMs", extensionBaMMs );
	} else {
		extensionBaMMs = strdup( "hbcp" );
	}

	if( opt >> GetOpt::OptionPresent( 'E', "extensionSeqs" ) ){
		opt >> GetOpt::Option( 'E', "extensionSeqs", extensionSeqs );
	} else {
		extensionSeqs = strdup( "fasta" );
	}

	if( opt >> GetOpt::OptionPresent( 'o', "outputDirectory" ) ){
		opt >> GetOpt::Option( 'o', "outputDirectory", outputDirectory );
		struct stat sb;
		if( stat( outputDirectory, &sb ) != 0 || !( S_ISDIR( sb.st_mode ) ) ){
			std::cerr << "Error: Output directory is not a directory" << std::endl;
		}
	}

	opt >> GetOpt::OptionPresent( 'v', "verbose", verbose );

	std::vector<std::string> argv;
	opt >> GetOpt::GlobalOption( argv );
	if( !( argv.size() == 2 ) ){
		std::cerr << "Error: Input directories are missing" << std::endl;
		printHelp();
		exit( -1 );
	}
	struct stat sb;
	inputDirectoryBaMMs = strdup( argv[0].c_str() );
	if( stat( inputDirectoryBaMMs, &sb ) != 0 || !( S_ISDIR( sb.st_mode ) ) ){
		std::cerr << "Error: Input directory with BaMM files is not a directory" << std::endl;
	}
	inputDirectorySeqs = strdup( argv[1].c_str() );
	if( stat( inputDirectorySeqs, &sb ) != 0 || !( S_ISDIR( sb.st_mode ) ) ){
		std::cerr << "Error: Input directory with FASTA files is not a directory" << std::endl;
	}

	if( outputDirectory == NULL ){
		outputDirectory = strdup( inputDirectorySeqs );
	}

	if( opt.options_remain() ){
		std::cerr << "Warning: Unknown arguments ignored" << std::endl;
	}

	return 0;
}

void Global::printHelp(){

	bool virusHostHelp = true; // print virus-host interaction or rather unspecific help comments
	bool developerHelp = false; // print developer-specific options

	printf( "\n" );
	printf( "SYNOPSIS\n" );
	if( virusHostHelp ){
		printf( "      predict BAMMDIRPATH VIRDIRPATH [OPTIONS]\n\n" );
	} else{
		printf( "      predict BAMMDIRPATH FASTADIRPATH [OPTIONS]\n\n" );
	}
	printf( "DESCRIPTION\n" );
	if( virusHostHelp ){
		printf( "      Predict virus-host interactions using pre-built homogeneous Bayesian\n"
				"      Markov models (BaMMs) of bacterial genomes.\n\n" );
		printf( "      BAMMDIRPATH\n"
				"          Input directory with homogeneous Bayesian Markov models (BaMMs) of\n"
				"          bacterial genomes. The default filename extension searched for is\n"
				"          <hbcp>. Use option -e (--extensionBaMMs) to change the default\n"
				"          setting.\n\n" );
		printf( "      VIRDIRPATH\n"
				"          Input directory with viral sequences in FASTA format (one file per\n"
				"          virus). The default filename extension searched for is <fasta>. Use\n"
				"          option -E (--extensionSeqs) to change the default setting.\n\n" );
	} else{
		printf( "      Calculate posterior probabilities for sequence sets in FASTA format using\n"
				"      pre-built homogeneous Bayesian Markov models (BaMMs).\n\n" );
		printf( "      BAMMDIRPATH\n"
				"          Input directory with homogeneous Bayesian Markov models (BaMMs). The\n"
				"          default filename extension searched for is <hbcp>. Use option -e\n"
				"          (--extensionBaMMs) to change the default setting.\n\n" );
		printf( "      FASTADIRPATH\n"
				"          Input directory with sequence sets in FASTA format. The default\n"
				"          filename extension searched for is <fasta>. Use option -E\n"
				"          (--extensionSeqs) to change the default setting.\n\n" );
	}
	printf( "OPTIONS\n" );
	printf( "  Options for input files\n" );
	printf( "      -e, --extensionBaMMs <STRING>\n"
			"          The filename extension of BaMM files searched for in BAMMDIRPATH. The\n"
			"          default is <hbcp>.\n\n" );
	if( virusHostHelp ){
		printf( "      -E, --extensionSeqs <STRING>\n"
				"          The filename extension of FASTA files searched for in VIRDIRPATH. The\n"
				"          default is <fasta>.\n\n" );
	} else{
		printf( "      -E, --extensionSeqs <STRING>\n"
				"          The filename extension of FASTA files searched for in FASTADIRPATH.\n"
				"          The default is <fasta>.\n\n" );
	}
	printf( "  Output options\n" );
	if( virusHostHelp ){
		printf( "      -o, --outputDirectory <DIRPATH>\n"
				"          Output directory for virus-host interactions. VIRDIRPATH is the\n"
				"          default output directory.\n\n" );
	} else{
		printf( "      -o, --outputDirectory <DIRPATH>\n"
				"          Output directory for posterior probabilities. FASTADIRPATH is the\n"
				"          default output directory.\n\n" );
	}
	printf( "      -v, --verbose\n"
			"          Verbose terminal printouts.\n\n" );
	printf( "      -h, --help\n"
			"          Printout this help.\n\n" );
	printf( "      --version\n"
			"          Show program name/version banner.\n\n" );
	if( developerHelp ){
		printf( "      (*) Developer options\n\n" );
	}
}
