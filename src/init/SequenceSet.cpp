#include "SequenceSet.h"
#include <boost/algorithm/string.hpp>

SequenceSet::SequenceSet( std::string sequenceFilepath,
                          bool singleStrand ){

	if( Alphabet::getSize() == 0 ){
		std::cerr << "Error: Initialize Alphabet before "
				"constructing a sequenceSet" << std::endl;
		exit( 1 );
	}

	sequenceFilepath_ = sequenceFilepath;

    baseSum_ = 0;

	baseFrequencies_ = new float[Alphabet::getSize()];

    isSingleStranded_ = singleStrand;

	readFASTA();
}

std::vector<std::string> SequenceSet::sequence2string(){
    std::vector<std::string> sequences;
    for( size_t n = 0; n< sequences_.size(); n++){
        std::string sequence;
        for(size_t i = 0; i < sequences_[n]->getL(); i++){
            sequence += Alphabet::getBase(sequences_[n]->getSequence()[i]);
        }
        sequences.push_back(sequence);
    }
    return sequences;
}

SequenceSet::~SequenceSet(){
    for( size_t i = 0; i < sequences_.size(); i++ ){
		delete sequences_[i];
	}
    
    delete[] baseFrequencies_;
}

std::string SequenceSet::getSequenceFilepath(){
	return sequenceFilepath_;
}

std::vector<Sequence*> SequenceSet::getSequences(){
	return sequences_;
}

size_t SequenceSet::getMinL(){
	return minL_;
}

size_t SequenceSet::getMaxL(){
	return maxL_;
}

size_t SequenceSet::getBaseSum(){
    return baseSum_;
}

float* SequenceSet::getBaseFrequencies(){
	return baseFrequencies_;
}

void SequenceSet::print(){

}

int SequenceSet::readFASTA(){

	/**
	 * while reading in the sequences do:
	 * 1. extract the header for each sequence
	 * 2. calculate the min. and max. sequence length in the file
	 * 3. count the number of sequences
	 * 4. calculate base frequencies
	 */

	size_t maxL = 0;
	size_t minL = std::numeric_limits<size_t>::max();
	std::vector<size_t> baseCounts( Alphabet::getSize() );
	std::string line, header, sequence;
	std::ifstream file( sequenceFilepath_.c_str() ); // opens FASTA file

	if( file.is_open() ){

		while( getline( file, line ) ){

			if( !( line.empty() ) ){ // skip blank lines

				if( line[0] == '>' ){

					if( !( header.empty() ) ){

						if( !( sequence.empty() ) ){

							size_t L = sequence.length();

                            maxL = ( L > maxL ) ? L : maxL;
                            minL = ( L < minL ) ? L : minL;

							// translate sequence into encoding
							uint8_t* encoding = ( uint8_t* )calloc( L, sizeof( uint8_t ) );

                            size_t n_count = 0;
							for( size_t i = 0; i < L; i++ ){

                                encoding[i] = Alphabet::getCode( sequence[i] );
								if( encoding[i] == 0 ){ // there is an unknown letter N
                                    n_count ++;
									continue;           // exclude undefined base from base counts
								}
                                baseCounts[encoding[i]-1]++; // count base
							}

                            if( n_count < L ) {
                                sequences_.push_back(new Sequence(encoding, L, header, isSingleStranded_));
                            } else{
                                std::cerr << "Warning: " << header
                                          << " contains only unknown alphabet and is filtered." << std::endl;
                            }
							sequence.clear();
							header.clear();

							free( encoding );

						} else {
							std::cerr << "Warning: Ignore FASTA entry without sequence: "
                                      << sequenceFilepath_ << std::endl;
							header.clear();
						}
					}

					if( line.length() == 1 ){ // corresponds to ">\n"
						// set header to sequence counter
						header = '>';
					} else {
						// fetch header till the first tab
                        // remove the '\r' terminator in the end of the line
                        std::vector<std::string> strs;
                        header = boost::split(strs, line, boost::is_any_of("\t"))[0];
                        header = boost::split(strs, header, boost::is_any_of("\r"))[0];
                    }

				} else if( !( header.empty() ) ){

					if( line.find( ' ' ) != std::string::npos ){
						// space character in sequence
						std::cerr << "Error: FASTA sequence contains space character: "
								  << sequenceFilepath_ << std::endl;
						exit( 1 );
					} else {
						sequence += line;
					}

				} else {
					std::cerr << "Error: Wrong FASTA format: " << sequenceFilepath_ << std::endl;
					exit( 1 );
				}
			}
		}

		if( !( header.empty() ) ){				// store the last sequence

			if( !( sequence.empty() ) ){

                size_t L = sequence.length();

                maxL = ( L > maxL ) ? L : maxL;
                minL = ( L < minL ) ? L : minL;

                // translate sequence into encoding
                uint8_t* encoding = ( uint8_t* )calloc( L, sizeof( uint8_t ) );

                size_t n_count = 0;
                for( size_t i = 0; i < L; i++ ){

                    encoding[i] = Alphabet::getCode( sequence[i] );

                    if( encoding[i] == 0 ){ // there is an unknown letter N
                        n_count++;
                        continue;           // exclude undefined base from base counts

                    }
                    baseCounts[encoding[i]-1]++; // count base
                }

                if( n_count < L ) {
                    sequences_.push_back(new Sequence(encoding, L, header, isSingleStranded_));
                } else{
                    std::cerr << "Warning: " << header
                              << " contains only unknown alphabet and is filtered." << std::endl;
                }

				sequence.clear();
				header.clear();

				free( encoding );

			} else {

				std::cerr << "Warning: Ignore FASTA entry without sequence: "
                          << sequenceFilepath_ << std::endl;
				header.clear();
			}
		}

		file.close();

	} else {

		std::cerr << "Error: Cannot open FASTA file: "
                  << sequenceFilepath_ << std::endl;
		exit( 1 );
	}

	maxL_ = maxL;
	minL_ = minL;

	 // calculate the sum of bases
	baseSum_ = 0;
	for( size_t i = 0; i < Alphabet::getSize(); i++ ){
		baseSum_ += baseCounts[i];
	}
    baseSum_ = (baseSum_ == 0) ? 1: baseSum_;

	// calculate base frequencies
	for( size_t i = 0; i < Alphabet::getSize(); i++ ){
		baseFrequencies_[i] = static_cast<float>( baseCounts[i] ) / static_cast<float>( baseSum_ );
	}

	return 0;
}
