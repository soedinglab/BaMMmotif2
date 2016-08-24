#ifndef GPREDICTION_H_
#define GPREDICTION_H_

#include <string>

#include <string.h>		// e.g. strdup
#include <sys/stat.h>

#include "../getopt_pp/getopt_pp.h"

#include "../shared/Alphabet.h"
#include "../shared/SequenceSet.h"

class Global{

public:

	static char*		inputDirectoryBaMMs;	// input directory with BaMM files
	static char*		inputDirectorySeqs;		// input directory with FASTA files

	static char*		outputDirectory;		// output directory

	static char*		extensionBaMMs;			// extension of BaMM files, defaults to hbcp
	static char*		extensionSeqs;			// extension of FASTA files, defaults to fasta

	static std::string	name;					// program name
	static std::string	version;				// program version number

	static bool			verbose;				// verbose printouts, defaults to false

	static void		init( int nargs, char* args[] );
	static void		destruct();

	static char* 	String( const char *s );			// GetOpt library

private:

	static int	    readArguments( int nargs, char* args[] );
	static void	    printHelp();
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

#endif /* GPREDICTION_H_ */
