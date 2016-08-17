#include "VirusHostInteractions.h"

VirusHostInteractions::VirusHostInteractions( char* inputDirectoryBaMMs, char* extensionBaMMs ){

	bamms_ = new BackgroundModelSet( inputDirectoryBaMMs, extensionBaMMs );

	for( size_t i = 0; i < bamms_->getN(); i++ ){
		bammNames_.push_back( bamms_->getBackgroundModels()[i]->getName() );
	}
}

VirusHostInteractions::~VirusHostInteractions(){

}

void VirusHostInteractions::predict( char* inputDirectorySeqs, char* extensionSeqs ){

	DIR* dir;
	struct dirent* ent;

	struct stat sb;

	if( ( dir = opendir( inputDirectorySeqs ) ) != NULL ){

		while( ( ent = readdir( dir ) ) != NULL ){

			std::string filep = std::string( inputDirectorySeqs ) + '/' + std::string( ent->d_name );

			if( stat( filep.c_str(), &sb ) == 0 && S_ISREG( sb.st_mode ) ){

				char* extension = strrchr( ent->d_name, '.' );
				if( strcmp( extension+1, extensionSeqs ) == 0 ){

					SequenceSet* sequenceSet = new SequenceSet( filep );
					sequenceSetNames_.push_back( baseName( sequenceSet->getSequenceFilepath().c_str() ) );

					std::vector<double> llikelihoods( bamms_->getN() );
					// calculate log likelihoods
					for( unsigned int i = 0; i < bamms_->getN(); i++ ){
						llikelihoods[i] = score( *bamms_->getBackgroundModels()[i], *sequenceSet );
					}
					// calculate posterior probabilities
					posteriors_.push_back( calculatePosteriors( llikelihoods ) );
				}
			}
		}
		closedir( dir );

		// calculate aggregate statistics
		aggregatePosteriors();

	} else{
		perror( inputDirectorySeqs );
	}

}

void VirusHostInteractions::print(){

	std::vector<size_t> bammNameIndices = sortIndices( bammNames_ );
	std::vector<size_t> sequenceSetNameIndices = sortIndices( sequenceSetNames_ );

	int wmax = 0;

	for( size_t i = 0; i < bammNames_.size(); i++ ){
		int w = static_cast<int>( bammNames_[i].length() );
		wmax = w > wmax ? w : wmax;
	}

	std::cout << std::setfill( ' ' ) << std::setw( wmax+1 ) << ' ';
	std::cout << sequenceSetNames_[sequenceSetNameIndices[0]];

	for( size_t j = 1; j < sequenceSetNames_.size(); j++ ){
		std::cout << ' ' << sequenceSetNames_[sequenceSetNameIndices[j]];
	}
	std::cout << std::endl;

	for( size_t i = 0; i < posteriors_.size(); i++ ){

		std::cout << std::setw( wmax );
		std::cout << bammNames_[bammNameIndices[i]];

		for( size_t j = 0; j < posteriors_[bammNameIndices[i]].size(); j++ ){

			std::cout << std::setw( static_cast<int>( sequenceSetNames_[sequenceSetNameIndices[j]].length()+1 ) );
			std::cout << std::fixed << std::setprecision( 3 );
			std::cout << posteriors_[bammNameIndices[i]][sequenceSetNameIndices[j]];
		}
		std::cout << std::endl;
	}
}

void VirusHostInteractions::write( char* outputDirectory ){

	std::ofstream file( std::string( outputDirectory ) + '/' + "vhi.tab" );
	if( file.is_open() ){

		std::vector<size_t> bammNameIndices = sortIndices( bammNames_ );
		std::vector<size_t> sequenceSetNameIndices = sortIndices( sequenceSetNames_ );

		file << "#";
		for( size_t j = 0; j < bestBammIndices_.size(); j++ ){
			file << " " << bammNameIndices[bestBammIndices_[sequenceSetNameIndices[j]]]+1;
		}
		file << std::endl;
		file << "#";
		for( size_t j = 0; j < bestBammPosteriors_.size(); j++ ){
			file << " " << std::fixed << std::setprecision( 2 ) << bestBammPosteriors_[sequenceSetNameIndices[j]];
		}
		file << std::endl;
		file << "#";
		for( size_t j = 0; j < secondBestBammPosteriors_.size(); j++ ){
			file << " " << std::fixed << std::setprecision( 2 ) << secondBestBammPosteriors_[sequenceSetNameIndices[j]];
		}
		file << std::endl;

		file << sequenceSetNames_[sequenceSetNameIndices[0]];
		for( size_t j = 1; j < sequenceSetNames_.size(); j++ ){
			file << '\t' << sequenceSetNames_[sequenceSetNameIndices[j]];
		}
		file << std::endl;

		for( size_t i = 0; i < posteriors_.size(); i++ ){

			file << bammNames_[bammNameIndices[i]];

			for( size_t j = 0; j < posteriors_[bammNameIndices[i]].size(); j++ ){

				file << '\t' << std::fixed << std::setprecision( 6 ) << posteriors_[bammNameIndices[i]][sequenceSetNameIndices[j]];
			}
			file << std::endl;
		}
		file.close();

	} else{

		std::cerr << "Error: Cannot write into output directory: " << outputDirectory << std::endl;
		exit( -1 );
	}
}

void VirusHostInteractions::aggregatePosteriors(){

	bestBammIndices_.clear();
	bestBammIndices_.resize( sequenceSetNames_.size(), 0 );

	bestBammPosteriors_.clear();
	bestBammPosteriors_.resize( sequenceSetNames_.size(), 0.0 );

	secondBestBammPosteriors_.clear();
	secondBestBammPosteriors_.resize( sequenceSetNames_.size(), 0.0 );

	for( size_t i = 0; i < bammNames_.size(); i++ ){
		for( size_t j = 0; j < sequenceSetNames_.size(); j++ ){

			if( posteriors_[i][j] > bestBammPosteriors_[j] ){

				secondBestBammPosteriors_[j] = bestBammPosteriors_[j];

				bestBammIndices_[j] = i;
				bestBammPosteriors_[j] = posteriors_[i][j];

			} else if( posteriors_[i][j] > secondBestBammPosteriors_[j] ){

				secondBestBammPosteriors_[j] = posteriors_[i][j];
			}
		}
	}
}

std::vector<double> VirusHostInteractions::calculatePosteriors( std::vector<double> llikelihoods ){

	// see http://stats.stackexchange.com/questions/66616/converting-normalizing-very-small-likelihood-values-to-probability/66621#66621

	int d = std::numeric_limits<double>::digits10; // digits of precision

	double epsilon = pow( 10, -d );
	long unsigned int N = llikelihoods.size();

	double limit = log( epsilon ) - log( static_cast<double>( N ) );

	// sort indices into ascending order
	std::vector<size_t> order = sortIndices( llikelihoods );
	// sort likelihoods into ascending order
	std::sort( llikelihoods.begin(), llikelihoods.end() );

	for( unsigned int i = 0; i < N; i++ ){
		llikelihoods[i] = llikelihoods[i] - llikelihoods[N-1];
	}

	double partition = 0.0;

	for( unsigned int i = 0; i < N; i++ ){
		if( llikelihoods[i] >= limit ){
			llikelihoods[i] = exp( llikelihoods[i] );
			partition += llikelihoods[i];
		}
	}

	std::vector<double> posteriors( N );

	for( unsigned int i = 0; i < N; i++ ){
		if( llikelihoods[i] >= limit ){
			posteriors[order[i]] = llikelihoods[i] / partition;
		} else{
			posteriors[order[i]] = 0.0;
		}
	}

	return( posteriors );
}

double VirusHostInteractions::score( BackgroundModel& bamm, SequenceSet& sequenceSet ){

	if( !( bamm.vIsLog() ) ){
		bamm.logV();
	}

	double likelihood = 0.0;

	int N = sequenceSet.getN();
	for( int n = 0; n < N; n++ ){

		int L = sequenceSet.getSequences()[n]->getL();
		for( int i = 0; i < L; i++ ){

			int k = std::min( i, bamm.getOrder() );
			int y = sequenceSet.getSequences()[n]->extractKmer( i, k );

			if( y >= 0 ){ // skip non-defined alphabet letters
				likelihood += bamm.getV()[k][y]; // add log probabilities
			}
		}
	}

	return likelihood;
}
