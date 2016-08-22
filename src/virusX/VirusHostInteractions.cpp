#include "VirusHostInteractions.h"

VirusHostInteractions::VirusHostInteractions( char* inputDirectoryBaMMs, char* extensionBaMMs ){

	bamms_ = new BackgroundModelSet( inputDirectoryBaMMs, extensionBaMMs );

	for( size_t i = 0; i < bamms_->getN(); i++ ){
		bammNames_.push_back( bamms_->getBackgroundModels()[i]->getName() );
	}
}

VirusHostInteractions::~VirusHostInteractions(){
	if( bamms_ != NULL ){
		delete bamms_;
	}
}

void VirusHostInteractions::predict( char* inputDirectorySeqs, char* extensionSeqs ){

	DIR* dir;
	struct dirent* ent;

	struct stat sb;

	if( ( dir = opendir( inputDirectorySeqs ) ) != NULL ){

		posteriors_.clear();
		posteriors_.resize( bamms_->getN() );

		while( ( ent = readdir( dir ) ) != NULL ){

			std::string filep = std::string( inputDirectorySeqs ) + '/' + std::string( ent->d_name );

			if( stat( filep.c_str(), &sb ) == 0 && S_ISREG( sb.st_mode ) ){

				char* extension = strrchr( ent->d_name, '.' );
				if( strcmp( extension+1, extensionSeqs ) == 0 ){

					SequenceSet* sequenceSet = new SequenceSet( filep );
					sequenceSetNames_.push_back( baseName( sequenceSet->getSequenceFilepath().c_str() ) );

//					// calculate log likelihoods
//					std::vector<double> llikelihoods( bamms_->getN() );
//					for( size_t i = 0; i < bamms_->getN(); i++ ){
//						llikelihoods[i] = score( *bamms_->getBackgroundModels()[i], *sequenceSet );
//					}
//
//					// calculate posterior probabilities
//					std::vector<double> posteriors = calculatePosteriors( llikelihoods );
//					for( size_t i = 0; i < posteriors_.size(); i++ ){
//						posteriors_[i].push_back( posteriors[i] );
//					}

					// calculate posterior probabilities
					std::vector<double> posteriors = bamms_->calculatePosteriorProbabilities( *sequenceSet );
					for( size_t i = 0; i < posteriors_.size(); i++ ){
						posteriors_[i].push_back( posteriors[i] );
					}
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

	std::cout << " _________" << std::endl;
	std::cout << "|*       *|" << std::endl;
	std::cout << "| Summary |" << std::endl;
	std::cout << "|*_______*|" << std::endl;
	std::cout << std::endl;

	int w1 = 0;
	for( size_t i = 0; i < sequenceSetNames_.size(); i++ ){
		int w = static_cast<int>( sequenceSetNames_[i].length() );
		w1 = w > w1 ? w : w1;
	}

	int w2 = 0;
	for( size_t i = 0; i < bestBammIndices_.size(); i++ ){
		int w = static_cast<int>( bammNames_[bestBammIndices_[i]].length() );
		w2 = w > w2 ? w : w2;
	}

	int w3 = 0;
	for( size_t i = 0; i < secondBestBammIndices_.size(); i++ ){
		int w = static_cast<int>( bammNames_[secondBestBammIndices_[i]].length() );
		w3 = w > w3 ? w : w3;
	}

	for( size_t j = 0; j < sequenceSetNames_.size(); j++ ){

		std::cout << std::setw( w1 );
		std::cout << sequenceSetNames_[sequenceSetNameIndices[j]];
		std::cout << '\t';
		std::cout << std::setw( w2+1 );
		std::cout << bammNames_[bestBammIndices_[sequenceSetNameIndices[j]]];
		std::cout << std::setw( 5 ) << std::fixed << std::setprecision( 2 );
		std::cout << bestBammPosteriors_[sequenceSetNameIndices[j]];
		std::cout << '\t';
		std::cout << std::setw( w3+1 );
		std::cout << bammNames_[secondBestBammIndices_[sequenceSetNameIndices[j]]];
		std::cout << std::setw( 5 ) << std::fixed << std::setprecision( 2 );
		std::cout << secondBestBammPosteriors_[sequenceSetNameIndices[j]];
		std::cout << std::endl;
	}

	std::cout << " _________________________" << std::endl;
	std::cout << "|*                       *|" << std::endl;
	std::cout << "| Virus-Host Interactions |" << std::endl;
	std::cout << "|*_______________________*|" << std::endl;
	std::cout << std::endl;

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
			std::cout << std::fixed << std::setprecision( 2 );
			std::cout << posteriors_[bammNameIndices[i]][sequenceSetNameIndices[j]];
		}
		std::cout << std::endl;
	}
}

void VirusHostInteractions::write( char* outputDirectory ){

	std::vector<size_t> bammNameIndices = sortIndices( bammNames_ );
	std::vector<size_t> sequenceSetNameIndices = sortIndices( sequenceSetNames_ );

	std::ofstream file( std::string( outputDirectory ) + '/' + "vhi.summary" );
	if( file.is_open() ){

		for( size_t j = 0; j < sequenceSetNames_.size(); j++ ){

			file << sequenceSetNames_[sequenceSetNameIndices[j]];
			file << '\t';
			file << bammNames_[bestBammIndices_[sequenceSetNameIndices[j]]];
			file << '\t';
			file << std::fixed << std::setprecision( 2 ) << bestBammPosteriors_[sequenceSetNameIndices[j]];
			file << '\t';
			file << bammNames_[secondBestBammIndices_[sequenceSetNameIndices[j]]];
			file << '\t';
			file << std::fixed << std::setprecision( 2 ) << secondBestBammPosteriors_[sequenceSetNameIndices[j]];
			file << std::endl;
		}
		file.close();

	} else{

		std::cerr << "Error: Cannot write into output directory: " << outputDirectory << std::endl;
		exit( -1 );
	}

	file.clear();
	file.open( std::string( outputDirectory ) + '/' + "vhi.matrix" );
	if( file.is_open() ){

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

	secondBestBammIndices_.clear();
	secondBestBammIndices_.resize( sequenceSetNames_.size(), 0 );

	secondBestBammPosteriors_.clear();
	secondBestBammPosteriors_.resize( sequenceSetNames_.size(), 0.0 );

	for( size_t i = 0; i < bammNames_.size(); i++ ){
		for( size_t j = 0; j < sequenceSetNames_.size(); j++ ){

			if( posteriors_[i][j] > bestBammPosteriors_[j] ){

				secondBestBammIndices_[j] = bestBammIndices_[j];
				secondBestBammPosteriors_[j] = bestBammPosteriors_[j];

				bestBammIndices_[j] = i;
				bestBammPosteriors_[j] = posteriors_[i][j];

			} else if( posteriors_[i][j] > secondBestBammPosteriors_[j] ){

				secondBestBammIndices_[j] = i;
				secondBestBammPosteriors_[j] = posteriors_[i][j];
			}
		}
	}
}
