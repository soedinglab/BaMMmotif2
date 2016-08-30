#ifndef UTILS_H_
#define UTILS_H_

#include <algorithm>	// e.g. std::sort
#include <limits>		// e.g. std::numeric_limits
#include <numeric>		// e.g. std::numeric
#include <vector>

#include <sys/stat.h>	// e.g. stat

static char*							baseName( const char* filePath );
// calculate posterior probabilities from log likelihoods
std::vector<double>						calculatePosteriorProbabilities( std::vector<double> lLikelihoods );
static void								createDirectory( char* dir );
static std::vector< std::vector<int> >	generateFoldIndices( unsigned int N, unsigned int folds );
// calculate the power for integer base
static int								ipow( unsigned int base, int exp );

template <typename T>
std::vector<size_t> sortIndices( const std::vector<T> &v ); // returns a permutation which rearranges v into ascending order

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

inline std::vector<double> calculatePosteriorProbabilities( std::vector<double> lLikelihoods ){

	// see http://stats.stackexchange.com/questions/66616/converting-normalizing-very-small-likelihood-values-to-probability/66621#66621

	int d = std::numeric_limits<double>::digits10; // digits of precision

	double epsilon = pow( 10, -d );
	long unsigned int N = lLikelihoods.size();

	double limit = log( epsilon ) - log( static_cast<double>( N ) );

	// sort indices into ascending order
	std::vector<size_t> order = sortIndices( lLikelihoods );
	// sort likelihoods into ascending order
	std::sort( lLikelihoods.begin(), lLikelihoods.end() );

	for( unsigned int i = 0; i < N; i++ ){
		lLikelihoods[i] = lLikelihoods[i] - lLikelihoods[N-1];
	}

	double partition = 0.0;

	for( unsigned int i = 0; i < N; i++ ){
		if( lLikelihoods[i] >= limit ){
			lLikelihoods[i] = exp( lLikelihoods[i] );
			partition += lLikelihoods[i];
		}
	}

	std::vector<double> posteriors( N );

	for( unsigned int i = 0; i < N; i++ ){
		if( lLikelihoods[i] >= limit ){
			posteriors[order[i]] = lLikelihoods[i] / partition;
		} else{
			posteriors[order[i]] = 0.0;
		}
	}

	return( posteriors );
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

template <typename T>
inline std::vector<size_t> sortIndices( const std::vector<T>& v ){

  // initialize with original indices
  std::vector<size_t> idx( v.size() );
  iota( idx.begin(), idx.end(), 0 );

  // sort indices based on comparing values in v
  sort( idx.begin(), idx.end(),
		[&v]( size_t i1, size_t i2 ){
	        return v[i1] < v[i2];
        }
      );

  return idx;
}

#endif /* UTILS_H_ */
