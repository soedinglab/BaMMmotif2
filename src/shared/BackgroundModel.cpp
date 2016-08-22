#include "BackgroundModel.h"

BackgroundModel::BackgroundModel( SequenceSet& sequenceSet,
		                          int order,
		                          std::vector<float> alpha,
		                          bool interpolate,
		                          std::vector< std::vector<int> > foldIndices,
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
	interpolate_ = interpolate;

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

	n_ = ( int** )malloc( ( K_+1 ) * sizeof( int* ) );
	for( int k = 0; k <= K_; k++ ){
		n_[k] = ( int* )calloc( Y_[k+1], sizeof( int ) );
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
			for( int k = 0; k <= K_; k++ ){
				// loop over sequence positions
				for( int i = k; i < L; i++ ){
					// extract (k+1)mer
					int y = sequenceSet.getSequences()[s_idx]->extractKmer( i, k );
					// skip non-defined alphabet letters
					if( y >= 0 ){
						// count (k+1)mer
						n_[k][y]++;

					}
				}
			}
		}
	}

	v_ = ( float** )malloc( ( K_+1 ) * sizeof( float* ) );
	for( int k = 0; k <= K_; k++ ){
		v_[k] = ( float* )calloc( Y_[k+1], sizeof( float ) );
	}
	// calculate conditional probabilities from counts
	calculateV();
}

BackgroundModel::BackgroundModel( std::string filePath ){

	name_.assign( baseName( filePath.c_str() ) );

	n_ = NULL;
	v_ = NULL;

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

			v_ = ( float** )malloc( ( K_+1 ) * sizeof( float* ) );
			for( int k = 0; k <= K_; k++ ){
				v_[k] = ( float* )calloc( Y_[k+1], sizeof( float ) );
			}

			float value;
			for( int k = 0; k <= K_; k++ ){
				for( int y = 0; y < Y_[k+1]; y++ ){

					if( fscanf( file, "%e", &value ) != EOF ){
						v_[k][y] = value;

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

	if( n_ != NULL ){
		for( int k = 0; k <= K_; k++ ){
			free( n_[k] );
		}
		free( n_ );
	}
	if( v_ != NULL ){
		for( int k = 0; k <= K_; k++ ){
			free( v_[k] );
		}
		free( v_ );
	}
}

std::string BackgroundModel::getName(){
    return name_;
}

int BackgroundModel::getOrder(){
    return K_;
}

float** BackgroundModel::getV(){
    return v_;
}

void BackgroundModel::expV(){

	for( int k = 0; k <= K_; k++ ){
		for( int y = 0; y < Y_[k+1]; y++ ){
			v_[k][y] = expf( v_[k][y] );
		}
	}
	vIsLog_ = false;
}

void BackgroundModel::logV(){

	for( int k = 0; k <= K_; k++ ){
		for( int y = 0; y < Y_[k+1]; y++ ){
			v_[k][y] = logf( v_[k][y] );
		}
	}
	vIsLog_ = true;
}

bool BackgroundModel::vIsLog(){
	return vIsLog_;
}

double BackgroundModel::calculateLogLikelihood( SequenceSet& sequenceSet,
		                                        std::vector<std::vector<int>> foldIndices,
		                                        std::vector<int> folds ){

	if( !( vIsLog_ ) ){
		logV();
	}

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

	double lLikelihood = 0.0;

	// loop over folds
	for( size_t f = 0; f < folds.size(); f++ ){
		// loop over fold indices
		for( size_t f_idx = 0; f_idx < foldIndices[folds[f]].size(); f_idx++ ){
			// get sequence index
			int s_idx = foldIndices[folds[f]][f_idx];
			// get sequence length
			int L = sequenceSet.getSequences()[s_idx]->getL();
			// loop over sequence positions
			for( int i = 0; i < L; i++ ){
				// calculate k
				int k = std::min( i, K_ );
				// extract (k+1)mer
				int y = sequenceSet.getSequences()[s_idx]->extractKmer( i, k );
				// skip non-defined alphabet letters
				if( y >= 0 ){
					// add log probabilities
					lLikelihood += v_[k][y];
				}
			}
		}
	}

	return lLikelihood;
}

void BackgroundModel::calculatePosLogLikelihoods( SequenceSet& sequenceSet,
		                                          char* outputDirectory,
		                                          std::vector<std::vector<int>> foldIndices,
		                                          std::vector<int> folds ){

	if( !( vIsLog_ ) ){
		logV();
	}

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

	std::ofstream file( std::string( outputDirectory ) + '/'
			            + name_ + '_'
			            + baseName( sequenceSet.getSequenceFilepath().c_str() )
			            + ".lhs" );

	if( file.is_open() ){

		// loop over folds
		for( size_t f = 0; f < folds.size(); f++ ){
			// loop over fold indices
			for( size_t f_idx = 0; f_idx < foldIndices[folds[f]].size(); f_idx++ ){
				// get sequence index
				int s_idx = foldIndices[folds[f]][f_idx];
				// get sequence length
				int L = sequenceSet.getSequences()[s_idx]->getL();
				// loop over sequence positions
				for( int i = 0; i < ( L-K_ ); i++ ){
					// reset log likelihood
					double lLikelihood = 0.0;
					// loop over order
					for( int k = 0; k <= K_; k++ ){
						// extract (k+1)mer
						int y = sequenceSet.getSequences()[s_idx]->extractKmer( i+k, k );
						if( y >= 0 ){
							// sum up log probabilities
							lLikelihood += v_[k][y];
						} else{
							// cannot calculate log probability for (k+1)mer and neither (K_+1)mer
							lLikelihood = 1.0;
							break;
						}
					}
					file << ( i == 0 ? "" : " " );
					if( lLikelihood > 0.0 ){
						file << "NA";
					} else{
						file << std::scientific << std::setprecision( 6 ) << lLikelihood;
					}
				}
				file << std::endl;
			}
		}
		file.close();

	} else{

		std::cerr << "Error: Cannot write into output directory: " << outputDirectory << std::endl;
		exit( -1 );
	}
}

void BackgroundModel::print(){

	if( interpolate_ ){
		std::cout << " ___________________________________" << std::endl;
		std::cout << "|*                                 *|" << std::endl;
		std::cout << "| Homogeneous Bayesian Markov Model |" << std::endl;
		std::cout << "|*_________________________________*|" << std::endl;
		std::cout << std::endl;
	} else{
		std::cout << " __________________________" << std::endl;
		std::cout << "|*                        *|" << std::endl;
		std::cout << "| Homogeneous Markov Model |" << std::endl;
		std::cout << "|*________________________*|" << std::endl;
		std::cout << std::endl;
	}

	std::cout << "name = " << name_ << std::endl << std::endl;
	std::cout << "K = " << K_ << std::endl;
	std::cout << "A =";
	for( int k = 0; k <= K_; k++ ){
		std::cout << " " << A_[k];
	}
	std::cout << std::endl;

	if( n_ != NULL ){

		std::cout << " ________" << std::endl;
		std::cout << "|        |" << std::endl;
		std::cout << "| Counts |" << std::endl;
		std::cout << "|________|" << std::endl;
		std::cout << std::endl;

		for( int k = 0; k <= K_; k++ ){
			std::cout << n_[k][0];
			for( int y = 1; y < Y_[k+1]; y++ ){
				std::cout << " " << n_[k][y];
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
		std::cout << std::fixed << std::setprecision( 3 ) << v_[k][0];
		for( int y = 1; y < Y_[k+1]; y++ ){
			std::cout << " " << std::fixed << std::setprecision( 3 ) << v_[k][y];
		}
		std::cout << std::endl;
	}
}

void BackgroundModel::write( char* dir ){

	std::ofstream file( std::string( dir ) + '/' + name_ + ( interpolate_ ? ".hb" : ".hnb" ) );
	if( file.is_open() ){

		file << "# K = " << K_ << std::endl;
		file << "# A =";
		for( int k = 0; k <= K_; k++ ){
			file << " " << A_[k];
		}
		file << std::endl;

		for( int k = 0; k <= K_; k++ ){
			file << std::scientific << std::setprecision( 6 ) << v_[k][0];
			for( int y = 1; y < Y_[k+1]; y++ ){
				file << std::scientific << std::setprecision( 6 ) << " " << v_[k][y];
			}
			file << std::endl;
		}
		file.close();

	} else{

		std::cerr << "Error: Cannot write into output directory: " << dir << std::endl;
		exit( -1 );
	}
}

void BackgroundModel::calculateV(){

	int baseCounts = 0;
	for( int y = 0; y < Y_[1]; y++ ){
		baseCounts += n_[0][y];
	}

	// calculate probabilities for order k = 0
	for( int y = 0; y < Y_[1]; y++ ){
		v_[0][y] = ( static_cast<float>( n_[0][y] ) + A_[0] * 0.25f ) // cope with absent bases using pseudocounts
				   / ( static_cast<float>( baseCounts ) + A_[0] );
	}

	// calculate probabilities for order k > 0
	for( int k = 1; k <= K_; k++ ){
		for( int y = 0; y < Y_[k+1]; y++ ){
			// omit first base (e.g. CGT) from y (e.g. ACGT)
			int y2 = y % Y_[k];
			// omit last base (e.g. ACG) from y (e.g. ACGT)
			int yk = y / Y_[1];
			if( interpolate_ ){
				v_[k][y] = ( static_cast<float>( n_[k][y] ) + A_[k] * v_[k-1][y2] )
						   / ( static_cast<float>( n_[k-1][yk] ) + A_[k] );
			} else{
				v_[k][y] = ( static_cast<float>( n_[k][y] ) + A_[k] * 0.25f )
						   / ( static_cast<float>( n_[k-1][yk] ) + A_[k] );
			}
		}
		// normalize probabilities
		float factor = 0.0f;
		for( int y = 0; y < Y_[k+1]; y++ ){
			factor += v_[k][y];
			if( ( y+1 ) % Y_[1] == 0 ){
				for( int i = 0; i < Y_[1]; i++ ){
					v_[k][y-i] /= factor;
				}
				factor = 0.0f;
			}
		}
	}
}
