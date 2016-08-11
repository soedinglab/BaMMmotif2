#ifndef UTILS_H_
#define UTILS_H_

#include <vector>

#include <sys/stat.h>	// e.g. stat

static char*							baseName( const char* filePath );
static void								createDirectory( char* dir );
static std::vector<std::vector<int>>	generateFolds( unsigned int N, unsigned int folds );
static int								ipow( unsigned int base, int exp );	// calculate the power for integer base

inline char* baseName( const char* filePath ){

	int i = 0, start = 0, end = 0;

	while( filePath[++i] != '\0' ){
		if( filePath[i] == '.' ){
			end = i - 1;
		}
	}
	while( --i != 0 && filePath[i] != '/' ){
		;
	}
	if( i != 0 ){
		start = i + 1;
	}

	char* baseName = ( char* )malloc( ( end-start+2 ) * sizeof( char ) );
	for( i = start; i <= end; i++ ){
		baseName[i-start] = filePath[i];
	}
	baseName[i-start] = '\0';

	return baseName;
}

inline void createDirectory( char* dir ){

	struct stat fileStatus;

	if( stat( dir, &fileStatus ) != 0 ){

		char* cmd = ( char* )calloc( 1024, sizeof( char ) );
		if( system( cmd ) != 0 ){
			fprintf( stderr, "Error: Directory %s could not be created\n", dir );
			exit( -1 );
		}
		free( cmd );
	}
}

inline std::vector<std::vector<int>> generateFolds( unsigned int N, unsigned int folds ){

	std::vector<std::vector<int>> indices( folds );

	for( unsigned int i = 0; i < N; i++ ){
		for( unsigned int j = 0; j < folds; j++ ){
			indices[j].push_back( i );
		}
	}

	return indices;
}

inline int ipow( unsigned int base, int exp ){

    int res = 1;

    while( exp ){
        if( exp & 1 )
        	res *= base;
        exp >>= 1;
        base *= base;
    }

    return res;
}



#endif /* UTILS_H_ */
