#include "BackgroundModel.h"

BackgroundModel::BackgroundModel( SequenceSet& sequenceSet,
		                          int order,
		                          std::vector<float> alpha,
		                          std::vector<std::vector<int>> foldIndices,
		                          std::vector<int> folds ){


	name_.assign( baseName( sequenceSet.getSequenceFilepath().c_str() ) );

	// calculate the maximum possible order
	int l = static_cast<int>( floorf(
				logf( static_cast<float>( std::numeric_limits<int>::max() ) ) /
				logf( static_cast<float>( Alphabet::getSize() ) ) ) );
	if( ( order+1 ) > l ){
		std::cerr << "Error: The maximum possible order is " << l << std::endl;
		exit( -1 );
	}

	K_ = order;
	for( int k = 0; k < K_+2; k++ ){
		Y_.push_back( ipow( Alphabet::getSize(), k ) );
	}

	A_ = alpha;

	if( folds.empty() ){
		if( foldIndices.empty() ){

			folds.push_back( 0 );

			foldIndices.push_back( std::vector<int>() );
			for( int n = 0; n < sequenceSet.getN(); n++ ){
				foldIndices[0].push_back( n );
			}
		} else{
			folds.resize( foldIndices.size() );
			std::iota( std::begin( folds ), std::end( folds ), 0 );
		}
	} else if( foldIndices.empty() ){

		folds.clear();
		folds.push_back( 0 );

		foldIndices.push_back( std::vector<int>() );
		for( int n = 0; n < sequenceSet.getN(); n++ ){
			foldIndices[0].push_back( n );
		}
	}

	n_bg_ = ( int** )malloc( ( K_+1 ) * sizeof( int* ) );
	for( int k = 0; k <= K_; k++ ){
		n_bg_[k] = ( int* )calloc( Y_[k+1], sizeof( int ) );
	}
	// calculate counts
	// loop over folds
	for( unsigned int f = 0; f < folds.size(); f++ ){
		// loop over fold indices
		for( unsigned int f_idx = 0; f_idx < foldIndices[folds[f]].size(); f_idx++ ){
			// get sequence index
			int s_idx = foldIndices[folds[f]][f_idx];
			// get sequence length
			int L = sequenceSet.getSequences()[s_idx]->getL();
			// loop over order
			for( int k = 0; k < ( K_+1 ); k++ ){
				// loop over sequence positions
				for( int i = k; i < L; i++ ){
					// extract (k+1)mer
					int y = sequenceSet.getSequences()[s_idx]->extractKmer( i, std::min( i, k ) );
					// skip non-defined alphabet letters
					if( y >= 0 ){
						// count (k+1)mer
						n_bg_[k][y]++;
					}
				}
			}
		}
	}

	v_bg_ = ( float** )malloc( ( K_+1 ) * sizeof( float* ) );
	for( int k = 0; k <= K_; k++ ){
		v_bg_[k] = ( float* )calloc( Y_[k+1], sizeof( float ) );
	}
	// calculate conditional probabilities from counts
	calculateVbg();
}

BackgroundModel::BackgroundModel( std::string filePath ){

	name_.assign( baseName( filePath.c_str() ) );

	n_bg_ = NULL;
	v_bg_ = NULL;

	struct stat sb;

	if( stat( filePath.c_str(), &sb ) == 0 && S_ISREG( sb.st_mode ) ){

		FILE* file;
		if( ( file = fopen( filePath.c_str(), "r" ) ) != NULL ){

			int K;
			if( fscanf( file, "# K = %d\n", &K ) == 1 ){
				K_ = K;
			} else{
				std::cerr << "Error: Wrong BaMM format: " << filePath << std::endl;
				exit( -1 );
			}
			A_.resize( K_+1 );

			float A;
			if( fscanf( file, "# A = %e", &A ) == 1 ){
				A_[0] = A;
			} else{
				std::cerr << "Error: Wrong BaMM format: " << filePath << std::endl;
				exit( -1 );
			}
			for( int k = 1; k <= K_; k++ ){
				if( fscanf( file, "%e", &A ) == 1 ){
					A_[k] = A;
				} else{
					std::cerr << "Error: Wrong BaMM format: " << filePath << std::endl;
					exit( -1 );
				}
			}

			for( int k = 0; k < K_+2; k++ ){
				Y_.push_back( ipow( Alphabet::getSize(), k ) );
			}

			v_bg_ = ( float** )malloc( ( K_+1 ) * sizeof( float* ) );
			for( int k = 0; k <= K_; k++ ){
				v_bg_[k] = ( float* )calloc( Y_[k+1], sizeof( float ) );
			}

			float value;
			for( int k = 0; k <= K_; k++ ){
				for( int y = 0; y < Y_[k+1]; y++ ){

					if( fscanf( file, "%e", &value ) != EOF ){
						v_bg_[k][y] = value;

					} else{

						std::cerr << "Error: Wrong BaMM format: " << filePath << std::endl;
						exit( -1 );
					}
				}
			}
			fclose( file );

		} else{

			std::cerr << "Error: Cannot open BaMM file: " << filePath << std::endl;
			exit( -1 );
		}
	}
}

BackgroundModel::~BackgroundModel(){

	if( n_bg_ != NULL ){
		for( int k = 0; k <= K_; k++ ){
			free( n_bg_[k] );
		}
		free( n_bg_ );
	}
	if( v_bg_ != NULL ){
		for( int k = 0; k <= K_; k++ ){
			free( v_bg_[k] );
		}
		free( v_bg_ );
	}
}

std::string BackgroundModel::getName(){
    return name_;
}

float** BackgroundModel::getVbg(){
    return v_bg_;
}

void BackgroundModel::print(){

	std::cout << " ___________________" << std::endl;
	std::cout << "|*                 *|" << std::endl;
	std::cout << "| Homogeneous BaMM! |" << std::endl;
	std::cout << "|*_________________*|" << std::endl;
	std::cout << std::endl;

	std::cout << "name = " << name_ << std::endl << std::endl;
	std::cout << "K = " << K_ << std::endl;
	std::cout << "A =";
	for( int k = 0; k <= K_; k++ ){
		std::cout << " " << A_[k];
	}
	std::cout << std::endl;

	if( n_bg_ != NULL ){

		std::cout << " ________" << std::endl;
		std::cout << "|        |" << std::endl;
		std::cout << "| Counts |" << std::endl;
		std::cout << "|________|" << std::endl;
		std::cout << std::endl;

		for( int k = 0; k <= K_; k++ ){
			std::cout << n_bg_[k][0];
			for( int y = 1; y < Y_[k+1]; y++ ){
				std::cout << " " << n_bg_[k][y];
			}
			std::cout << std::endl;
		}
	}

	std::cout << " ___________________________" << std::endl;
	std::cout << "|                           |" << std::endl;
	std::cout << "| Conditional probabilities |" << std::endl;
	std::cout << "|___________________________|" << std::endl;
	std::cout << std::endl;

	for( int k = 0; k <= K_; k++ ){
		std::cout << std::fixed << std::setprecision( 3 ) << v_bg_[k][0];
		for( int y = 1; y < Y_[k+1]; y++ ){
			std::cout << " " << std::fixed << std::setprecision( 3 ) << v_bg_[k][y];
		}
		std::cout << std::endl;
	}
}

void BackgroundModel::write( char* dir ){

	std::ofstream file( std::string( dir ) + '/' + name_ + ".hb" );

	file << "# K = " << K_ << std::endl;
	file << "# A =";
	for( int k = 0; k <= K_; k++ ){
		file << " " << A_[k];
	}
	file << std::endl;

	for( int k = 0; k <= K_; k++ ){
		file << std::fixed << std::setprecision( 6 ) << v_bg_[k][0];
		for( int y = 1; y < Y_[k+1]; y++ ){
			file << std::fixed << std::setprecision( 6 ) << " " << v_bg_[k][y];
		}
		file << std::endl;
	}
}

void BackgroundModel::calculateVbg(){

	int baseCounts = 0;
	for( int y = 0; y < Y_[1]; y++ ){
		baseCounts += n_bg_[0][y];
	}

	// calculate probabilities for order k = 0
	for( int y = 0; y < Y_[1]; y++ ){
		v_bg_[0][y] = ( static_cast<float>( n_bg_[0][y] ) + A_[0] * 0.25f ) / // cope with absent bases using pseudocounts
				      ( static_cast<float>( baseCounts ) + A_[0] );
	}

	// calculate probabilities for order k > 0
	for( int k = 1; k <= K_; k++ ){
		for( int y = 0; y < Y_[k+1]; y++ ){
			// omit first base (e.g. CGT) from y (e.g. ACGT)
			int y2 = y % Y_[k];
			// omit last base (e.g. ACG) from y (e.g. ACGT)
			int yk = y / Y_[1];
			v_bg_[k][y] = ( static_cast<float>( n_bg_[k][y] ) + A_[k] * v_bg_[k-1][y2] ) /
					      ( static_cast<float>( n_bg_[k-1][yk] ) + A_[k] );
		}
		// normalize probabilities
		float factor = 0.0f;
		for( int y = 0; y < Y_[k+1]; y++ ){
			factor += v_bg_[k][y];
			if( ( y+1 ) % Y_[1] == 0 ){
				for( int i = 0; i < Y_[1]; i++ ){
					v_bg_[k][y-i] /= factor;
				}
				factor = 0.0f;
			}
		}
	}
}
