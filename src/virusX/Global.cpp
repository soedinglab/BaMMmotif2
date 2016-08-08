#include "Global.h"

char*               Global::inputDirectory = NULL;				// input directory
char*               Global::outputDirectory = NULL;				// output directory

char*				Global::extension = NULL;					// extension of files in FASTA format, defaults to fasta

int					Global::modelOrder = 2;						// background model order
std::vector<float>	Global::modelAlpha( modelOrder+1, 1.0f );	// background model alpha
float				Global::modelBeta = 10.0f;					// alpha_k = beta x gamma^(k-1) for k > 0
float				Global::modelGamma = 2.0f;					// - " -

bool                Global::verbose = false;					// verbose printouts


void Global::init( int nargs, char* args[] ){

	readArguments( nargs, args );

	char* alphabetString = strdup( "STANDARD" );
	Alphabet::init( alphabetString );
	free( alphabetString );
}

void Global::destruct(){

	if( inputDirectory != NULL ){
		free( inputDirectory );
	}
	if( outputDirectory != NULL ){
		free( outputDirectory );
	}
	if( extension != NULL ){
		free( extension );
	}

	Alphabet::destruct();
}

int Global::readArguments( int nargs, char* args[] ){

	GetOpt::GetOpt_pp opt( nargs, args );

	if( opt >> GetOpt::OptionPresent( 'h', "help" ) ){
		printHelp();
		exit( -1 );
	}

	if( opt >> GetOpt::OptionPresent( 'e', "extension" ) ){
		opt >> GetOpt::Option( 'e', "extension", extension );
	} else {
		extension = strdup( "fasta" );
	}

	opt >> GetOpt::Option( 'k', modelOrder );

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
		float b;
		if( opt >> GetOpt::OptionPresent( 'b', "beta" ) ){
			opt >> GetOpt::Option( 'b', "beta", b );
			fprintf( stderr, "Option -b is ignored.\n" );
		}
		float g;
		if( opt >> GetOpt::OptionPresent( 'g', "gamma" ) ){
			opt >> GetOpt::Option( 'g', "gamma", g );
			fprintf( stderr, "Option -g is ignored.\n" );
		}
	} else{
		if( opt >> GetOpt::OptionPresent( 'b', "beta" ) ){
			opt >> GetOpt::Option( 'b', "beta", modelBeta );
		}
		if( opt >> GetOpt::OptionPresent( 'g', "gamma" ) ){
			opt >> GetOpt::Option( 'g', "gamma", modelGamma );
		}
		if( static_cast<int>( modelAlpha.size() ) != modelOrder+1 ){
			if( static_cast<int>( modelAlpha.size() ) > modelOrder+1 ){
				modelAlpha.resize( modelOrder+1 );
			} else{
				modelAlpha.resize( modelOrder+1, modelAlpha.back() );
			}
		}
		if( modelOrder > 0 ){
			for( unsigned int i=1; i < modelAlpha.size(); i++ ){
				modelAlpha[i] = modelBeta * powf( modelGamma, static_cast<float>( i )-1.0f );
			}
		}
	}

	if( opt >> GetOpt::OptionPresent( 'o', "outputDirectory" ) ){
		opt >> GetOpt::Option( 'o', "outputDirectory", outputDirectory );
	}

	opt >> GetOpt::OptionPresent( 'v', "verbose", verbose );

	std::vector<std::string> argv;
	opt >> GetOpt::GlobalOption( argv );
	if( !( argv.size() == 1 ) ){
		std::cerr << "Error: Path to directory with genomes in FASTA format missing" << std::endl;
		printHelp();
		exit( -1 );
	}
	inputDirectory = strdup( argv[0].c_str() );

	if( outputDirectory == NULL ){
		outputDirectory = strdup( inputDirectory );
	}

	if( opt.options_remain() ){
		std::cerr << "Warning: Unknown arguments ignored" << std::endl;
	}

	return 0;
}

void Global::printHelp(){

	bool developerHelp = false; // print developer-specific options

	printf( "\n" );
	printf( "SYNOPSIS\n" );
	printf( "      VirusX DIRPATH [OPTIONS]\n\n" );
	printf( "DESCRIPTION\n" );
	printf( "      Learn Bayesian homogeneous Markov models from bacterial genomes.\n\n" );
	printf( "      DIRPATH\n"
			"          Input directory with bacterial genomes in FASTA format (one file per\n"
			"          genome). The default filename extension searched for is .fasta. Use\n"
			"          option -e (--extension) to change the default setting.\n\n" );
	printf("OPTIONS\n");
	printf( "  Options for input files\n" );
	printf( "      -e, --extension <STRING>\n"
			"          The filename extension of FASTA files searched for in the input\n"
			"          directory. The default is <fasta>.\n\n" );
	printf( "  Options for Bayesian homogeneous Markov models\n" );
	printf( "      -k, --order <INTEGER>\n"
			"          Order. The default is 2.\n\n" );
	printf( "      -a, --alpha <FLOAT> [<FLOAT>...] (*)\n"
			"          Order-specific prior strength. The default is 1.0 (for k = 0) and\n"
			"          20 x 3^(k-1) (for k > 0). The options -b and -g are ignored.\n\n" );
	printf( "      -b, --beta <FLOAT> (*)\n"
			"          Calculate order-specific alphas according to beta x gamma^(k-1) (for\n"
			"          k > 0). The default is 20.0.\n\n" );
	printf( "      -g, --gamma <FLOAT> (*)\n"
			"          Calculate order-specific alphas according to beta x gamma^(k-1) (for\n"
			"          k > 0). The default is 3.0.\n\n" );
	printf( "  Output options\n" );
	printf( "      -o, --outputDirectory <DIRPATH>\n"
			"          Output directory for Bayesian homogenous Markov models. The input\n"
			"          directory is the default output directory.\n\n" );
	printf( "      -v, --verbose\n"
			"          Verbose terminal printouts.\n\n" );
	printf( "      -h, --help\n"
			"          Printout this help.\n\n" );
	if( developerHelp ){
		printf( "      (*) Developer options\n\n" );
	}
}
