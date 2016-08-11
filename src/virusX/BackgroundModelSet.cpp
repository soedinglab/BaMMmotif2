#include "BackgroundModelSet.h"

BackgroundModelSet::BackgroundModelSet( char* inputDirectory, char* extension, int order, std::vector<float> alpha ){

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
					BackgroundModel* bamm = new BackgroundModel( *sequenceSet, order, alpha );

					backgroundModels_.push_back( bamm );
					N_++;
				}
			}
		}
		closedir( dir );

	} else{
		perror( Global::inputDirectory );
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
					N_++;
				}
			}
		}
		closedir( dir );

	} else{
		perror( Global::inputDirectory );
	}
}

BackgroundModelSet::~BackgroundModelSet(){

	for( std::list<BackgroundModel*>::iterator iter = backgroundModels_.begin(); iter != backgroundModels_.end(); ){
		delete *iter;
		iter = backgroundModels_.erase( iter );
	}
}

std::list<BackgroundModel*>& BackgroundModelSet::getBackgroundModels(){
    return backgroundModels_;
}

int BackgroundModelSet::getN(){
    return N_;
}

void BackgroundModelSet::print(){

	for( std::list<BackgroundModel*>::iterator iter = backgroundModels_.begin(); iter != backgroundModels_.end(); iter++ ){
		( *iter )->print();
	}
}

void BackgroundModelSet::write( char* dir ){

	for( std::list<BackgroundModel*>::iterator iter = backgroundModels_.begin(); iter != backgroundModels_.end(); iter++ ){
		( *iter )->write( dir );
	}
}
