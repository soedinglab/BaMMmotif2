#include <fstream>		// std::fstream

#include "MotifSet.h"

MotifSet::MotifSet(){

	if( Global::BaMMpatternFilename != NULL ){

		// scan file and conduct for IUPAC pattern
		std::ifstream file;
		file.open( Global::BaMMpatternFilename, std::ifstream::in );
		std::string pattern;
		int length;
		if( !file.good() ){
			std::cout << "Error: Cannot open pattern sequence file: "
					<< Global::BaMMpatternFilename << std::endl;
			exit( -1 );
		} else {
			while( getline( file, pattern ) ){

				length = ( int ) pattern.length();			// get length of the first sequence

				length += Global::addColumns.at( 0 ) + Global::addColumns.at( 1 );

				Motif* motif = new Motif( length );

				motif->initFromBaMMPattern( pattern );

				motifs_.push_back( motif );
			}
		}

	} else if( Global::bindingSiteFilename != NULL ){

		std::ifstream file;
		file.open( Global::bindingSiteFilename, std::ifstream::in );
		std::string bindingSite;
		int length;

		if( !file.good() ){
			std::cout << "Error: Cannot open bindingSitesFile sequence file: "
					<< Global::bindingSiteFilename << std::endl;
			exit( -1 );
		} else {
			getline( file, bindingSite );					// get length of the first sequence
			length = ( int ) bindingSite.length();
		}

		length += Global::addColumns.at( 0 ) + Global::addColumns.at( 1 );

		Motif* motif = new Motif( length );

		motif->initFromBindingSites( Global::bindingSiteFilename );

		motifs_.push_back( motif );

		N_ = 1;

	} else if( Global::PWMFilename != NULL ){

		// read file to calculate motif length
		std::ifstream file;
		file.open( Global::PWMFilename, std::ifstream::in );

		if( !file.good() ){
			std::cout << "Error: Cannot open PWM file: " << Global::PWMFilename << std::endl;
			exit( -1 );
		} else {

			int length;									// length of motif with extension
			int asize;									// alphabet size
			std::string line;
			std::string row;

			int y, j;
			while( getline( file, line ) ){
				// search for the line starting with "letter-probability matrix"
				if( line[0] == 'l' ){

					// get the size of alphabet after "alength= "
					std::stringstream Asize( line.substr( line.find( "h=" ) + 2 ) );
					Asize >> asize;

					// get the length of motif after "w= "
					std::stringstream Width( line.substr( line.find( "w=" ) + 2 ) );
					Width >> length;

					// extend the length due to addColumns option
					length += Global::addColumns.at( 0 ) + Global::addColumns.at( 1 );

					// construct an initial motif
					Motif* motif = new Motif( length );

					// parse PWM for each motif
					float** PWM = new float*[asize];

					for( y = 0; y < asize; y++ ){

						PWM[y] = new float[length];

						for( j = 0; j < Global::addColumns.at( 0 ); j++ ){
							PWM[y][j] = 1.0f / ( float )asize;
						}
						for( j = length - Global::addColumns.at( 1 ); j < length; j++ ){
							PWM[y][j] = 1.0f / ( float )asize;
						}

					}

					// get the following W lines
					for( j = Global::addColumns.at( 0 ); j < length - Global::addColumns.at( 1 ) ; j++ ){

						getline( file, row );

						std::stringstream number( row );

						float weight;

						y = 0;
						while( number >> weight ){
							PWM[y][j] = weight;
							y++;
						}

					}

					// initialize each motif with a PWM
					motif->initFromPWM( PWM, asize, N_ );

					motifs_.push_back( motif );

					// count the number of motifs
					N_ ++;

					// free memory for PWM
					for( y = 0; y < asize; y++ ){
						delete PWM[y];
					}
					delete PWM;
				}
			}
		}
	} else if( Global::BaMMFilename != NULL ){

		// read file to calculate motif length
		std::ifstream file;
		file.open( Global::BaMMFilename, std::ifstream::in );

		if( !file.good() ){
			std::cout << "Error: Cannot open BaMM file: " << Global::BaMMFilename << std::endl;
			exit( -1 );
		} else {

			int length = 0;
			int order  = -1;
			int asize  = 0;
			std::string line;


			// read file to calculate motif length, order and alphabet size
			while( getline( file, line ) ){

				if( line.empty() ){
					length++;
				}else{
					if( length == 0 ){
						order++;
						if( order == 0 ){
							for( unsigned int i = 0; i < line.length(); i++ ){
								if( line[i] == ' ' )
									asize++;
							}
						}
					}
				}
			}

			// adjust modelOrder
			Global::modelOrder = order;

			// adjust Alphabet
			char* alphabetType = new char[9];
			if( asize == 4 ){
				strcpy( alphabetType, "STANDARD" );
			} else {
				if( asize == 6 ){
					strcpy( alphabetType, "EXTENDED" );
				}else{
					// asize isn't enough to guess alphabet type: throw error
					std::cout << "Error: Please provide Alphabet Type of your BaMM File: " << Global::BaMMFilename << std::endl;
					exit( -1 );
				}
			}
			// todo: what if asize = 5?

			Alphabet::init( alphabetType );

			// construct an initial motif
			Motif* motif = new Motif( length );

			// initialize motif from file
			motif->initFromBaMM( Global::BaMMFilename );

			motifs_.push_back( motif );

			N_ = 1;
		}
	}
	//	if( Global::verbose ) print();
}

MotifSet::~MotifSet(){
    for( size_t i = 0; i < motifs_.size(); i++ ){
        delete motifs_[i];
    }
}

std::vector<Motif*> MotifSet::getMotifs(){
    return motifs_;
}

int MotifSet::getN(){
    return N_;
}

void MotifSet::print(){
	printf( " ____________________________\n"
			"|                            |\n"
			"| PROBABILITIES for MotifSet |\n"
			"|____________________________|\n\n" );

	for( int i = 0; i < N_; i++ ){
		printf( " ________________________________________\n"
				"|                                        |\n"
				"| INITIALIZED PROBABILITIES for Motif %d  |\n"
				"|________________________________________|\n\n", i+1 );
		motifs_[i]->print();
	}
}
void MotifSet::write(){

}
