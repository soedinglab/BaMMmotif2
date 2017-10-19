#include "BackgroundModelSet.h"

BackgroundModelSet::BackgroundModelSet( char* indir, char* extension,
				size_t order, std::vector<float> alpha, bool interpolate ){

	DIR* dir;

	struct dirent* ent;

	struct stat sb;

	if( ( dir = opendir( indir ) ) != NULL ){

		while( ( ent = readdir( dir ) ) != NULL ){

			std::string filep = std::string( indir ) + '/'
								+ std::string( ent->d_name );

			if( stat( filep.c_str(), &sb ) == 0 && S_ISREG( sb.st_mode ) ){

				char* ext = strrchr( ent->d_name, '.' );
				if( strcmp( ext+1, extension ) == 0 ){

					SequenceSet* sequenceSet = new SequenceSet( filep );

					BackgroundModel* bamm = new BackgroundModel(
										sequenceSet->getSequences(),
										order, alpha, interpolate,
										baseName( sequenceSet->getSequenceFilepath().c_str() ) );

					backgroundModels_.push_back( bamm );
					// todo: here the memory for sequenceSet and bamm is not freed
				}
			}
		}
		closedir( dir );

	} else {
		perror( indir );
        exit( 1 );
	}
}

BackgroundModelSet::BackgroundModelSet( char* indir, char* extension ){

	DIR* dir;
	struct dirent* ent;

	struct stat sb;

	if( ( dir = opendir( indir ) ) != NULL ){

		while( ( ent = readdir( dir ) ) != NULL ){

			std::string filep = std::string( indir ) + '/' +
								std::string( ent->d_name );

			if( stat( filep.c_str(), &sb ) == 0 && S_ISREG( sb.st_mode ) ){

				char* ext = strrchr( ent->d_name, '.' );
				if( strcmp( ext+1, extension ) == 0 ){

					BackgroundModel* bamm = new BackgroundModel( filep );
					backgroundModels_.push_back( bamm );
				}
			}
		}
		closedir( dir );

	} else {
		perror( indir );
        exit( 1 );
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

std::vector<double> BackgroundModelSet::calculateLogLikelihoods(
		std::vector<Sequence*> seqs ){

	std::vector<double> llikelihoods( backgroundModels_.size() );

	for( size_t i = 0; i < backgroundModels_.size(); i++ ){

		llikelihoods[i] = backgroundModels_[i]->calculateLogLikelihood( seqs );
	}
	return llikelihoods;
}

std::vector<double> BackgroundModelSet::calculatePosteriorProbabilities(
		std::vector<Sequence*> seqs ){

	return ::calculatePosteriorProbabilities( calculateLogLikelihoods( seqs ) );
}

void BackgroundModelSet::calculatePosLikelihoods( std::vector<Sequence*> seqs,
		char* odir ){

	for( size_t i = 0; i < backgroundModels_.size(); i++ ){
		backgroundModels_[i]->calculatePosLikelihoods( seqs, odir );
	}
}

void BackgroundModelSet::print(){

/*	for( size_t i = 0; i < backgroundModels_.size(); i++ ){
		backgroundModels_[i]->print();
	}*/
}

void BackgroundModelSet::write( char* odir ){

	for( size_t i = 0; i < backgroundModels_.size(); i++ ){
		backgroundModels_[i]->write( odir, backgroundModels_[i]->getName() );
	}
}
