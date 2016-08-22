#include "BackgroundModelSet.h"

BackgroundModelSet::BackgroundModelSet( char* inputDirectory, char* extension, int order, std::vector<float> alpha, bool interpolate ){

	DIR* dir;
	struct dirent* ent;

	struct stat sb;

	if( ( dir = opendir( inputDirectory ) ) != NULL ){

		while( ( ent = readdir( dir ) ) != NULL ){

			std::string filep = std::string( inputDirectory ) + '/' + std::string( ent->d_name );

			if( stat( filep.c_str(), &sb ) == 0 && S_ISREG( sb.st_mode ) ){

				char* ext = strrchr( ent->d_name, '.' );
				if( strcmp( ext+1, extension ) == 0 ){

					SequenceSet* sequenceSet = new SequenceSet( filep );
					BackgroundModel* bamm = new BackgroundModel( *sequenceSet, order, alpha, interpolate );

					backgroundModels_.push_back( bamm );
				}
			}
		}
		closedir( dir );

	} else{
		perror( inputDirectory );
	}
}

BackgroundModelSet::BackgroundModelSet( char* inputDirectory, char* extension ){

	DIR* dir;
	struct dirent* ent;

	struct stat sb;

	if( ( dir = opendir( inputDirectory ) ) != NULL ){

		while( ( ent = readdir( dir ) ) != NULL ){

			std::string filep = std::string( inputDirectory ) + '/' + std::string( ent->d_name );

			if( stat( filep.c_str(), &sb ) == 0 && S_ISREG( sb.st_mode ) ){

				char* ext = strrchr( ent->d_name, '.' );
				if( strcmp( ext+1, extension ) == 0 ){

					BackgroundModel* bamm = new BackgroundModel( filep );
					backgroundModels_.push_back( bamm );
				}
			}
		}
		closedir( dir );

	} else{
		perror( inputDirectory );
	}
}

BackgroundModelSet::~BackgroundModelSet(){

	for( size_t i = 0; i < backgroundModels_.size(); i++ ){
		delete backgroundModels_[i];
	}
}

std::vector<BackgroundModel*>& BackgroundModelSet::getBackgroundModels(){
    return backgroundModels_;
}

size_t BackgroundModelSet::getN(){
    return( backgroundModels_.size() );
}

std::vector<double> BackgroundModelSet::calculateLogLikelihoods( SequenceSet& sequenceSet ){

	std::vector<double> llikelihoods( backgroundModels_.size() );

	for( size_t i = 0; i < backgroundModels_.size(); i++ ){

		llikelihoods[i] = backgroundModels_[i]->calculateLogLikelihood( sequenceSet );
	}
	return llikelihoods;
}

std::vector<double> BackgroundModelSet::calculatePosteriorProbabilities( SequenceSet& sequenceSet ){

	return ::calculatePosteriorProbabilities( calculateLogLikelihoods( sequenceSet ) );
}

void BackgroundModelSet::print(){

	for( size_t i = 0; i < backgroundModels_.size(); i++ ){
		backgroundModels_[i]->print();
	}
}

void BackgroundModelSet::write( char* dir ){

	for( size_t i = 0; i < backgroundModels_.size(); i++ ){
		backgroundModels_[i]->write( dir );
	}
}
