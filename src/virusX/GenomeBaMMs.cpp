#include "GenomeBaMMs.h"

GenomeBaMMs::GenomeBaMMs(){

	DIR* dir;
	struct dirent* ent;

	struct stat sb;

	if( ( dir = opendir( Global::inputDirectory ) ) != NULL ){

		while( ( ent = readdir( dir ) ) != NULL ){

			std::string filep = std::string( Global::inputDirectory ) + '/' + std::string( ent->d_name );

			if( stat( filep.c_str(), &sb ) == 0 && S_ISREG( sb.st_mode ) ){

				char* extension = strrchr( ent->d_name, '.' );
				if( strcmp( extension+1, Global::extension ) == 0 ){

					std::cout << ent->d_name << std::endl;
					SequenceSet* sequenceSet = new SequenceSet( filep );
					BackgroundModel* bamm = new BackgroundModel( *sequenceSet, Global::modelOrder, Global::modelAlpha );
					bamm->write( Global::outputDirectory, baseName( ent->d_name ) );
					bamms_.push_back( bamm );
				}
			}
		}
		closedir( dir );

	} else{
		perror( Global::inputDirectory );
	}
}

GenomeBaMMs::~GenomeBaMMs(){

}

void GenomeBaMMs::print(){

	for( std::list<BackgroundModel*>::iterator iter = bamms_.begin(); iter != bamms_.end(); iter++ ){
		( *iter )->print();
	}
}
