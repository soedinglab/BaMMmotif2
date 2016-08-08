#include "GenomeBaMMs.h"

GenomeBaMMs::GenomeBaMMs(){

	DIR* dir;
	struct dirent* ent;


	if( ( dir = opendir( Global::inputDirectory ) ) != NULL ){

		while( ( ent = readdir( dir ) ) != NULL ){

			char* extension = strrchr( ent->d_name, '.' );
			if( strcmp( extension+1, Global::extension ) == 0 ){

				std::string p = std::string( Global::inputDirectory ) + '/' + std::string( ent->d_name );
				std::cout << p << std::endl;
				SequenceSet* sequenceSet = new SequenceSet( p );
				sequenceSet->print();
//				BackgroundModel* bamm = new BackgroundModel( *sequenceSet, Global::modelOrder, Global::modelAlpha );
//				bamms_.push_back( bamm );
//				bamm->print();
			}
		}
		closedir( dir );
	} else {
		perror( Global::inputDirectory );
	}

}

GenomeBaMMs::~GenomeBaMMs(){

}
