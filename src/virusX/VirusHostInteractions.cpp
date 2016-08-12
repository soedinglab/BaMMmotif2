#include "VirusHostInteractions.h"

VirusHostInteractions::VirusHostInteractions( char* inputDirectoryBaMMs, char* extensionBaMMs ){

	bamms_ = new BackgroundModelSet( inputDirectoryBaMMs, extensionBaMMs );
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

					for( std::list<BackgroundModel*>::iterator iter = bamms_->getBackgroundModels().begin(); iter != bamms_->getBackgroundModels().end(); iter++ ){
						score( **iter, *sequenceSet );
					}

				}
			}
		}
		closedir( dir );

	} else{
		perror( inputDirectorySeqs );
	}

}

void VirusHostInteractions::print(){

}

void VirusHostInteractions::write( char* outputDirectory ){

}

float VirusHostInteractions::score( BackgroundModel& bamm, SequenceSet& sequenceSet ){
	std::cout << bamm.getName() << std::endl;
	std::cout << baseName( sequenceSet.getSequenceFilepath().c_str() ) << std::endl;
	return 0.0f;
}
