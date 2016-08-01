#ifndef GLOBAL_H_
#define GLOBAL_H_

#include <iostream>											// std::cout, std::endl
#include <vector>
#include <list>
#include <string>
#include <cstring> 											// std::strcpy
#include <numeric>											// std::iota
#include <iomanip>											// std::setprecision
#include <fstream>
#include <limits>											// std::numeric_limits<int>::max();

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <dirent.h> // directory entries
#include <sys/stat.h> // file status

#include "../getopt_pp/getopt_pp.h" // get options

#include "../shared/Alphabet.h"
#include "../shared/SequenceSet.h"

class Global{

public:

	static char*		inputDirectory; 					// input directory
	static char*		outputDirectory; 					// output directory

	static char*		negSequenceFilename;				// filename of negative sequence FASTA file

	static char*		negSequenceBasename;				// basename of negative sequence FASTA file

	static char*		extension;							// extensions of files in FASTA format, default is fasta

	static char* 		alphabetString;						// provide alphabet type
	static bool			revcomp;							// also search on reverse complement of sequences, defaults to false

	static SequenceSet*	negSequenceSet;						// negative Sequence Set

	// background model options
	static int			bgModelOrder;						// background model order, defaults to 2
	static float		bgModelAlpha;						// background model alpha

	static std::vector< std::vector<int> >	negFoldIndices;	// sequence indices for each cross-validation fold

	static bool			verbose;							// verbose printouts, defaults to false

	static void         init( int nargs, char* args[] );
	static void         destruct();

	static char* 		String( const char *s );			// convert const char* to string, for GetOpt library
	static int			ipow( unsigned int base, int exp );	// power function for integers
	static int*			powA;								// sizes of alphabet to the power k

private:

	static int	        readArguments( int nargs, char* args[] );
	static void	        createDirectory( char* dir );
	static char*		baseName( char* path );
	static void	        generateFolds();
	static void	        printHelp();
};

inline char* Global::String( const char *s ){
	return strdup( s );
}

namespace GetOpt{
	template <> inline _Option::Result convert<char*>( const std::string& s, char*& d, std::ios::fmtflags ){
		_Option::Result ret = _Option::BadType;
		d = Global::String( s.c_str() );
		ret = _Option::OK;
		return ret;
	}
}

#endif /* GLOBAL_H_ */
