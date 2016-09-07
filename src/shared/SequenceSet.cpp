#include "SequenceSet.h"

SequenceSet::SequenceSet( std::string sequenceFilepath, bool revcomp, std::string intensityFilepath ){

	if( Alphabet::getSize() == 0 ){
		std::cerr << "Error: Initialize Alphabet before constructing a SequenceSet" << std::endl;
		exit( -1 );
	}

	sequenceFilepath_ = sequenceFilepath;

	// calculate the length of the subsequence that is max. extractable
	// depends on the return type (int) of Sequence::extractKmer()
	int l = static_cast<int>( floorf(
				logf( static_cast<float>( std::numeric_limits<int>::max() ) ) /
				logf( static_cast<float>( Alphabet::getSize() ) ) ) );

	for( int i = 0; i <= l; i++ ){
		Y_.push_back( ipow( Alphabet::getSize(), i ) );
	}

	readFASTA( revcomp );

	if( !( intensityFilepath.empty() ) ){
		intensityFilepath_ = intensityFilepath;
		readIntensities();
	}
}

SequenceSet::~SequenceSet(){

    for( size_t i = 0; i < sequences_.size(); i++ ){
		delete sequences_[i];
	}

	if( baseFrequencies_ ){
		delete[] baseFrequencies_;
	}
}

std::string SequenceSet::getSequenceFilepath(){
	return sequenceFilepath_;
}

std::string SequenceSet::getIntensityFilepath(){
	return intensityFilepath_;
}

std::vector<Sequence*> SequenceSet::getSequences(){
	return sequences_;
}

int SequenceSet::getN(){
	return N_;
}

unsigned int SequenceSet::getMinL(){
	return minL_;
}

unsigned int SequenceSet::getMaxL(){
	return maxL_;
}

float* SequenceSet::getBaseFrequencies(){
	return baseFrequencies_;
}

void SequenceSet::print(){

	for( int i = 0; i < N_; i++ ){
		sequences_[i]->print();
	}
}

void SequenceSet::debug(){
    // read FASTA with massive printouts to check the function
    // -> NOTE: this is a "quick-an-dirty" version which is invalid as soon as someone changes something in readFATSA.
    // needs to be done better:
    //-> since readFASTA is independent of global variables, one should do a real test sequenceSet and check if
    // the output for readFasta gives the expected outcome instead

        bool revcomp = false;
        int N = 0; // sequence counter
        int maxL = 0;
        int minL = std::numeric_limits<int>::max();
        std::vector<unsigned int> baseCounts( Alphabet::getSize() );

        std::string line, header, sequence;
        std::ifstream file( sequenceFilepath_.c_str() ); // opens FASTA file

        if( file.is_open() ){

            while( getline( file, line ).good() ){

                fprintf( stdout,  "Current Line: %s \n", line.c_str());

                if( !( line.empty() ) ){ // skip blank lines

                    fprintf( stdout,  "Line not empty\n");

                    if( line[0] == '>' ){

                        fprintf( stdout,  "Line is a header line\n");

                        if( !( header.empty() ) ){

                            fprintf( stdout,  "There has already been a header line before\n");

                            if( !( sequence.empty() ) ){

                                fprintf( stdout,  "There has already been a sequence line before\n");

                                N++; // increment sequence counter
                                fprintf( stdout,   "\t N -> %d \n", N);

                                int L = static_cast<int>( sequence.length() );

                                fprintf( stdout,   "current sequence length -> %d \n", L);

                                if ( L > maxL ){
                                    fprintf( stdout, "maxL = %d \n ", maxL);
                                    fprintf( stdout,  "sequence length is longer than maxL -> maxL = L \n");
                                    maxL = L;
                                    fprintf( stdout, "maxL = %d \n ", maxL);
                                }
                                if ( L < minL ){
                                    fprintf( stdout, "minL = %d \n ", minL);
                                    fprintf( stdout,  "sequence length is smaller than minL -> minL = L \n");
                                    minL = L;
                                    fprintf( stdout, "minL = %d \n ", minL);
                                }

                                // translate sequence into encoding
                                uint8_t* encoding = ( uint8_t* )calloc( L, sizeof( uint8_t ) );
                                fprintf( stdout, "translate sequence into encoding ...\n");
                                fprintf( stdout, " Sequence ->  Encoding \n");
                                for( int i = 0; i < L; i++ ){
                                    encoding[i] = Alphabet::getCode( sequence[i] );
                                    fprintf( stdout, "       %d ->  %d       \n", sequence[i], encoding[i] );
                                    if( encoding[i] == 0 ){
                                        std::cerr << "Warning: The FASTA file contains an undefined base: " << sequence[i] << std::endl;
                                        continue; // exclude undefined base from base counts
                                    }
                                    baseCounts[encoding[i]-1]++; // count base
                                }

                                fprintf( stdout, "\n Counted baseCounts:\n");
                                fprintf( stdout, "Encoding -> counts \n");
                                for( int i = 0; i < L; i++ ){
                                    encoding[i] = Alphabet::getCode( sequence[i] );
                                    if( encoding[i] == 0 ){
                                        std::cerr << "Warning: The FASTA file contains an undefined base: " << sequence[i] << std::endl;
                                        continue; // exclude undefined base from base counts
                                    }
                                    fprintf( stdout, "      %d -> %d     \n", encoding[i], baseCounts[encoding[i]-1]);
                                }

                                fprintf( stdout, "create new Sequence with encoding, L, header, Y, revcomp");
                                sequences_.push_back( new Sequence( encoding, L, header, Y_, revcomp ) );

                                sequence.clear();
                                header.clear();

                                free( encoding );

                            } else{

                                std::cerr << "Warning: Ignore FASTA entry without sequence: " << sequenceFilepath_ << std::endl;
                                header.clear();
                            }
                        }

                        if( line.length() == 1 ){ // corresponds to ">\n"
                            // set header to sequence counter
                            header = static_cast<std::ostringstream*>( &( std::ostringstream() << ( N+1 ) ) )->str();
                        } else{
                            header = line.substr( 1 );// fetch header
                        }
                        fprintf( stdout, "Header = %s \n", header.c_str() );

                    } else if ( !( header.empty() ) ){

                        fprintf( stdout, "Line is a sequence and header is not empty\n");

                        if( line.find( ' ' ) != std::string::npos ){
                            // space character in sequence
                            std::cerr << "Error: FASTA sequence contains space character: " << sequenceFilepath_ << std::endl;
                            exit( -1 );
                        } else{
                            sequence += line;
                            fprintf( stdout, "line added to sequence\n");
                            fprintf( stdout, "%s", sequence.c_str() );
                        }

                    } else{

                        std::cerr << "Error: Wrong FASTA format: " << sequenceFilepath_ << std::endl;
                        exit( -1 );
                    }
                }
            }

            if( !( header.empty() ) ){ // store the last sequence

                if( !( sequence.empty() ) ){

                    N++; // increment sequence counter
                    fprintf( stdout,   "\t N -> %d \n", N);

                    int L = static_cast<int>( sequence.length() );
                    fprintf( stdout,   "current sequence length -> %d \n", L);

                    if ( L > maxL ){
                        fprintf( stdout, "maxL = %d \n ", maxL);
                        fprintf( stdout,  "sequence length is longer than maxL -> maxL = L \n");
                        maxL = L;
                        fprintf( stdout, "maxL = %d \n ", maxL);                    }
                    if ( L < minL ){
                        fprintf( stdout, "minL = %d \n ", minL);
                        fprintf( stdout,  "sequence length is smaller than minL -> minL = L \n");
                        minL = L;
                        fprintf( stdout, "minL = %d \n ", minL);
                    }

                    // translate sequence into encoding
                    uint8_t* encoding = ( uint8_t* )calloc( L, sizeof( uint8_t ) );
                    fprintf( stdout, "translate sequence into encoding ...\n");
                    fprintf( stdout, " Sequence ->  Encoding \n");
                    for( int i = 0; i < L; i++ ){
                        encoding[i] = Alphabet::getCode( sequence[i] );
                        fprintf( stdout, "       %d ->  %d       \n", sequence[i], encoding[i] );
                        if( encoding[i] == 0 ){
                            std::cerr << "Warning: The FASTA file contains an undefined base: " << sequence[i] << std::endl;
                            continue; // exclude undefined base from base counts
                        }
                        baseCounts[encoding[i]-1]++; // count base
                    }

                    fprintf( stdout, "\n Counted baseCounts:\n");
                    fprintf( stdout, "Encoding -> counts \n");
                    for( int i = 0; i < L; i++ ){
                        encoding[i] = Alphabet::getCode( sequence[i] );
                        if( encoding[i] == 0 ){
                            std::cerr << "Warning: The FASTA file contains an undefined base: " << sequence[i] << std::endl;
                            continue; // exclude undefined base from base counts
                        }
                        fprintf( stdout, "      %d -> %d     \n", encoding[i], baseCounts[encoding[i]-1]);
                    }

                    fprintf( stdout, "create new Sequence with encoding, L, header, Y, revcomp");
                    sequences_.push_back( new Sequence( encoding, L, header, Y_, revcomp ) );

                    sequence.clear();
                    header.clear();

                    free( encoding );

                } else{

                    std::cerr << "Warning: Ignore FASTA entry without sequence: " << sequenceFilepath_ << std::endl;
                    header.clear();
                }
            }

            file.close();

        } else{

            std::cerr << "Error: Cannot open FASTA file: " << sequenceFilepath_ << std::endl;
            exit( -1 );
        }

        N_ = N;
        maxL_ = maxL;
        minL_ = minL;

        fprintf( stdout, "Final Numbers: N_ = %d maxL_ = %d minL_ = %d \n", N, maxL, minL);

         // calculate the sum of bases
        unsigned int sumCounts = 0;
        for( int i = 0; i < Alphabet::getSize(); i++ )
            sumCounts += baseCounts[i];

        fprintf( stdout, "calculate the sum of bases = %d \n", sumCounts);
        fprintf( stdout, "calculate base frequencies =");
        // calculate base frequencies
        baseFrequencies_ = new float[Alphabet::getSize()];
        for( int i = 0; i < Alphabet::getSize(); i++ ){
            baseFrequencies_[i] = static_cast<float>( baseCounts[i] ) /
                    static_cast<float>( sumCounts );
            fprintf( stdout, " %.2f ", baseFrequencies_[i]);
        }
        fprintf( stdout, "\n\n");
        exit(0);

}

int SequenceSet::readFASTA( bool revcomp ){

	/**
	 * while reading in the sequences do
	 * 1. extract the header for each sequence
	 * 2. calculate the min. and max. sequence length in the file
	 * 3. count the number of sequences
	 * 4. calculate base frequencies
	 */

	int N = 0; // sequence counter
	int maxL = 0;
	int minL = std::numeric_limits<int>::max();
	std::vector<unsigned int> baseCounts( Alphabet::getSize() );

	std::string line, header, sequence;
	std::ifstream file( sequenceFilepath_.c_str() ); // opens FASTA file

	if( file.is_open() ){

		while( getline( file, line ).good() ){

			if( !( line.empty() ) ){ // skip blank lines

				if( line[0] == '>' ){

					if( !( header.empty() ) ){

						if( !( sequence.empty() ) ){

							N++; // increment sequence counter

							int L = static_cast<int>( sequence.length() );

							if ( L > maxL ){
								maxL = L;
							}
							if ( L < minL ){
								minL = L;
							}

							// translate sequence into encoding
							uint8_t* encoding = ( uint8_t* )calloc( L, sizeof( uint8_t ) );
							for( int i = 0; i < L; i++ ){
								encoding[i] = Alphabet::getCode( sequence[i] );
								if( encoding[i] == 0 ){
									std::cerr << "Warning: The FASTA file contains an undefined base: " << sequence[i] << std::endl;
									continue; // exclude undefined base from base counts
								}
								baseCounts[encoding[i]-1]++; // count base
							}
							sequences_.push_back( new Sequence( encoding, L, header, Y_, revcomp ) );

							sequence.clear();
							header.clear();

							free( encoding );

						} else{

							std::cerr << "Warning: Ignore FASTA entry without sequence: " << sequenceFilepath_ << std::endl;
							header.clear();
						}
					}

					if( line.length() == 1 ){ // corresponds to ">\n"
						// set header to sequence counter
						header = static_cast<std::ostringstream*>( &( std::ostringstream() << ( N+1 ) ) )->str();
					} else{
						header = line.substr( 1 );// fetch header
					}

				} else if ( !( header.empty() ) ){

					if( line.find( ' ' ) != std::string::npos ){
						// space character in sequence
						std::cerr << "Error: FASTA sequence contains space character: " << sequenceFilepath_ << std::endl;
						exit( -1 );
					} else{
						sequence += line;
					}

				} else{

					std::cerr << "Error: Wrong FASTA format: " << sequenceFilepath_ << std::endl;
					exit( -1 );
				}
			}
		}

		if( !( header.empty() ) ){ // store the last sequence

			if( !( sequence.empty() ) ){

				N++; // increment sequence counter

				int L = static_cast<int>( sequence.length() );

				if ( L > maxL ){
					maxL = L;
				}
				if ( L < minL ){
					minL = L;
				}

				// translate sequence into encoding
				uint8_t* encoding = ( uint8_t* )calloc( L, sizeof( uint8_t ) );
				for( int i = 0; i < L; i++ ){
					encoding[i] = Alphabet::getCode( sequence[i] );
					if( encoding[i] == 0 ){
						std::cerr << "Warning: The FASTA file contains an undefined base: " << sequence[i] << std::endl;
						continue; // exclude undefined base from base counts
					}
					baseCounts[encoding[i]-1]++; // count base
				}
				sequences_.push_back( new Sequence( encoding, L, header, Y_, revcomp ) );

				sequence.clear();
				header.clear();

				free( encoding );

			} else{

				std::cerr << "Warning: Ignore FASTA entry without sequence: " << sequenceFilepath_ << std::endl;
				header.clear();
			}
		}

		file.close();

	} else{

		std::cerr << "Error: Cannot open FASTA file: " << sequenceFilepath_ << std::endl;
		exit( -1 );
	}

	N_ = N;
	maxL_ = maxL;
	minL_ = minL;

	 // calculate the sum of bases
	unsigned int sumCounts = 0;
	for( int i = 0; i < Alphabet::getSize(); i++ )
		sumCounts += baseCounts[i];

	// calculate base frequencies
	baseFrequencies_ = new float[Alphabet::getSize()];
	for( int i = 0; i < Alphabet::getSize(); i++ )
		baseFrequencies_[i] = static_cast<float>( baseCounts[i] ) /
		                      static_cast<float>( sumCounts );

	return 0;
}

int SequenceSet::readIntensities(){

	std::cerr << "Error: SequenceSet::readIntensities() is not implemented so far." << std::endl;
	exit( -1 );
}
