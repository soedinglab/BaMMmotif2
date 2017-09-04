#include "MotifSet.h"

MotifSet::MotifSet( char* indir, size_t l_flank, size_t r_flank,
                    std::string tag, SequenceSet* posSet,
                    float* f_bg, size_t K, std::vector<float> alphas ){
	N_ = 0;

	if( tag.compare( "bindingsites" ) == 0 ){

		std::ifstream file;
		file.open( indir, std::ifstream::in );
		std::string bindingSite;
		size_t length;

		if( !file.good() ){
			std::cout << "Error: Cannot open binding sites file: "
					<< indir << std::endl;
			exit( -1 );
		} else {
			getline( file, bindingSite );	// get length of the first sequence
			length = bindingSite.length();
		}

		length += l_flank + r_flank;

		Motif* motif = new Motif( length, K, alphas, f_bg );

		motif->initFromBindingSites( indir, l_flank, r_flank );

		motifs_.push_back( motif );

		N_ = 1;

		// todo: here to delete motif
		// delete motif;  // Error

	} else if( tag.compare( "PWM" ) == 0 ){

		// read file to calculate motif length
		std::ifstream file;
		file.open( indir, std::ifstream::in );

		if( !file.good() ){

			std::cout << "Error: Cannot open PWM file: " << indir << std::endl;

			exit( -1 );

		} else {

			size_t length;						// length of motif
			size_t asize;						// alphabet size
			std::string line;
			std::string row;

			while( getline( file, line ) ){

				// search for the line starting with "letter-probability matrix"
				if( line[0] == 'l' ){

					// get the size of alphabet after "alength= "
					std::stringstream A( line.substr( line.find( "h=" ) + 2 ) );
					A >> asize;

					// get the length of motif after "w= "
					std::stringstream W( line.substr( line.find( "w=" ) + 2 ) );
					W >> length;

					// extend the length due to addColumns option
					length += l_flank + r_flank;

					// construct an initial motif
					Motif* motif = new Motif( length, K, alphas, f_bg );

					// parse PWM for each motif
					float** PWM = new float*[asize];

					for( size_t y = 0; y < asize; y++ ){

						PWM[y] = new float[length];

						for( size_t j = 0; j < l_flank; j++ ){
							PWM[y][j] = 1.0f / ( float )asize;
						}
						for( size_t j = length - r_flank; j < length; j++ ){
							PWM[y][j] = 1.0f / ( float )asize;
						}

					}

					// get the following W lines
					for( size_t j = l_flank; j < length - r_flank ; j++ ){

						getline( file, row );

						std::stringstream number( row );

						for( size_t y = 0; y < asize; y++ ){
							number >> PWM[y][j];
						}

					}

					// initialize each motif with a PWM
					motif->initFromPWM( PWM, asize, posSet );

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
		file.open( indir, std::ifstream::in );

		if( !file.good() ){

			std::cout << "Error: Cannot open file: " << indir << std::endl;
			exit( -1 );

		} else {

			size_t model_length = 0;
			size_t model_order = 0;
			std::string line;

			// read file to calculate motif length, order and alphabet size
			while( getline( file, line ) ){

				if( line.empty() ){
					// count the number of empty lines as the motif length
					model_length++;

				} else if( model_length == 0 ){
					// count the lines of the first position
					model_order++;

				}
			}

			// extend the core region of the model due to the added columns
			model_length += l_flank + r_flank;
			// adjust model order, extra 1
			model_order -= 1;

			// construct an initial motif
			Motif* motif = new Motif( model_length, K, alphas, f_bg );

			// initialize motif from file
			motif->initFromBaMM( indir, l_flank, r_flank );

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
