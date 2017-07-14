#include "MotifSet.h"

#include <fstream>		// std::fstream
#include <string>		// compare

MotifSet::MotifSet( char* indir, size_t l_flank, size_t r_flank,
		size_t order, std::string tag ){

	indir_ = indir;
	l_flank_ = l_flank;
	r_flank_ = r_flank;
	K_ = order;
	tag_ = tag;

	if( tag_.compare( "bindingsites" ) == 0 ){

		std::ifstream file;
		file.open( indir_, std::ifstream::in );
		std::string bindingSite;
		size_t length;

		if( !file.good() ){
			std::cout << "Error: Cannot open binding sites file: "
					<< indir_ << std::endl;
			exit( -1 );
		} else {
			getline( file, bindingSite );	// get length of the first sequence
			length = bindingSite.length();
		}

		length += l_flank_ + r_flank_;

		Motif* motif = new Motif( length );

		motif->initFromBindingSites( indir_ );

		motifs_.push_back( motif );

		N_ = 1;

	} else if( tag_.compare( "PWM" ) == 0 ){

		// read file to calculate motif length
		std::ifstream file;
		file.open( indir_, std::ifstream::in );
		size_t minL = Global::posSequenceSet->getMinL();

		if( !file.good() ){

			std::cout << "Error: Cannot open PWM file: " << indir_ << std::endl;

			exit( -1 );

		} else {

			size_t length;									// length of motif with extension
			size_t asize;									// alphabet size
			std::string line;
			std::string row;

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
					length += l_flank_ + r_flank_;

					if( length > minL ){
						fprintf( stderr, "Error: Width of PWM exceeds the length "
								"of the shortest sequence in the sequence set.\n" );
						exit( -1 );
					}

					// construct an initial motif
					Motif* motif = new Motif( length );

					// parse PWM for each motif
					float** PWM = new float*[asize];

					for( size_t y = 0; y < asize; y++ ){

						PWM[y] = new float[length];

						for( size_t j = 0; j < l_flank_; j++ ){
							PWM[y][j] = 1.0f / ( float )asize;
						}
						for( size_t j = length - r_flank_; j < length; j++ ){
							PWM[y][j] = 1.0f / ( float )asize;
						}

					}

					// get the following W lines
					for( size_t j = l_flank_; j < length - r_flank_ ; j++ ){

						getline( file, row );

						std::stringstream number( row );

						for( size_t y = 0; y < asize; y++ ){
							number >> PWM[y][j];
						}

					}

					// initialize each motif with a PWM
					motif->initFromPWM( PWM, asize );

					motifs_.push_back( motif );

					// count the number of motifs
					N_++;

					// free memory for PWM
					for( size_t y = 0; y < asize; y++ ){
						delete[] PWM[y];
					}
					delete[] PWM;
				}
			}
		}

	} else if( tag.compare( "BaMM" ) == 0 ){

		// each BaMM file contains one optimized motif model
		// read file to calculate motif length
		std::ifstream file;
		file.open( indir_, std::ifstream::in );

		if( !file.good() ){

			std::cout << "Error: Cannot open bamm file: " << indir_ << std::endl;
			exit( -1 );

		} else {

			size_t model_length = 0;
			size_t model_order_plus_1 = 0;
			std::string line;

			// read file to calculate motif length, order and alphabet size
			while( getline( file, line ) ){

				if( line.empty() ){
					// count the number of empty lines as the motif length
					model_length++;

				} else if( model_length == 0 ){
					// count the lines of the first position on the motif as the model model_order_plus_1
					model_order_plus_1++;

				}
			}

			// extend the core region of the model due to the added columns
			model_length += l_flank_ + r_flank_;

			// adjust model order due to the .bamm file
			Global::modelOrder = model_order_plus_1 - 1;

			// construct an initial motif
			Motif* motif = new Motif( model_length );

			// initialize motif from file
			motif->initFromBaMM( indir_ );

			motifs_.push_back( motif );

			N_ = 1;
		}
	}
}

MotifSet::~MotifSet(){
    for( size_t i = 0; i < motifs_.size(); i++ ){
        delete motifs_[i];
    }
}

std::vector<Motif*> MotifSet::getMotifs(){
    return motifs_;
}

size_t MotifSet::getN(){
    return N_;
}

void MotifSet::print(){

	fprintf(stderr, " ____________________________\n"
					"|                            |\n"
					"| PROBABILITIES for MotifSet |\n"
					"|____________________________|\n\n" );


	for( size_t i = 0; i < N_; i++ ){
		fprintf(stderr, " ________________________________________\n"
						"|                                        |\n"
						"| INITIALIZED PROBABILITIES for Motif %d  |\n"
						"|________________________________________|\n\n",
						( int )i+1 );
		motifs_[i]->print();
	}
}
void MotifSet::write( char* outdir ){

}
