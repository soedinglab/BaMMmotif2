#ifndef UTILS_H_
#define UTILS_H_

#include <vector>

#include <sys/stat.h>	// e.g. stat

static char*							baseName( const char* filePath );
static void								createDirectory( char* dir );
static int								ipow( unsigned int base, int exp );							// calculate the power for integer base
static std::vector< std::vector<int> >	generateFoldIndices( unsigned int N, unsigned int folds );
static void								quickSort( std::vector<float> arr, int left, int right );	// sort in descending order using

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

inline std::vector< std::vector<int> > generateFoldIndices( unsigned int N, unsigned int folds ){

	std::vector< std::vector<int> > indices( folds );

	for( unsigned int i = 0; i < N; i += folds ){
		for( unsigned int j = 0; j < folds; j++ ){
		    if( i+j < N ){
				indices[j].push_back( i+j );
		    }
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

inline void quickSort( std::vector<float> arr, int left, int right ){

	int i = left, j = right;
	float tmp;
	float pivot = arr[( left + right ) / 2];

	/* partition */
	while( i <= j ){
		while( arr[i] - pivot > 0 )	i++;
		while( arr[j] - pivot < 0 )	j--;
		if( i <= j ){
			tmp = arr[i];
			arr[i] = arr[j];
			arr[j] = tmp;
			i++;
			j--;
		}
	}

	/* recursion */
	if( left < j )	quickSort( arr, left, j );
	if( i < right )	quickSort( arr, i, right );
}

#endif /* UTILS_H_ */
