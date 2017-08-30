#include "SequenceSet.h"

SequenceSet::SequenceSet( std::string sequenceFilepath,
							bool singleStrand,
							std::string intensityFilepath,
							float q ){

	if( Alphabet::getSize() == 0 ){
		std::cerr << "Error: Initialize Alphabet before "
				"constructing a SequenceSet" << std::endl;
		exit( -1 );
	}

	sequenceFilepath_ = sequenceFilepath;

	for( size_t i = 0; i <= 11; i++ ){
		Y_.push_back( ipow( Alphabet::getSize(), i ) );
	}

	baseFrequencies_ = new float[Y_[1]];

	readFASTA( singleStrand );

	if( !( intensityFilepath.empty() ) ){
		intensityFilepath_ = intensityFilepath;
		readIntensities();
	}

	q_ = q;

}

SequenceSet::~SequenceSet(){
    for( size_t i = 0; i < sequences_.size(); i++ ){
		delete sequences_[i];
	}
    if( !baseFrequencies_ ) delete baseFrequencies_;
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

size_t SequenceSet::getMinL(){
	return minL_;
}

size_t SequenceSet::getMaxL(){
	return maxL_;
}

float SequenceSet::getQ(){
	return q_;
}

float* SequenceSet::getBaseFrequencies(){
	return baseFrequencies_;
}

void SequenceSet::print(){

}

int SequenceSet::readFASTA( bool singleStrand ){

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

		while( getline( file, line ).good() ){

			if( !( line.empty() ) ){ // skip blank lines

				if( line[0] == '>' ){

					if( !( header.empty() ) ){

						if( !( sequence.empty() ) ){

							size_t L = sequence.length();

							if( L > maxL ){
								maxL = L;
							}
							if( L < minL ){
								minL = L;
							}

							// translate sequence into encoding
							uint8_t* encoding = ( uint8_t* )calloc( L, sizeof( uint8_t ) );

							for( size_t i = 0; i < L; i++ ){

								encoding[i] = Alphabet::getCode( sequence[i] );

								if( encoding[i] == 0 ){

//									std::cerr << "Warning: The FASTA file contains an undefined base: " <<
//											sequence[i] << " at sequence " << header << std::endl;

									continue; // exclude undefined base from base counts

								}
								baseCounts[encoding[i]-1]++; // count base
							}
							sequences_.push_back( new Sequence( encoding, L, header, Y_, singleStrand ) );

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
						header = line.substr( 1 );  // fetch header
					}

				} else if( !( header.empty() ) ){

					if( line.find( ' ' ) != std::string::npos ){
						// space character in sequence
						std::cerr << "Error: FASTA sequence contains space character: "
								<< sequenceFilepath_ << std::endl;
						exit( -1 );
					} else {
						sequence += line;
					}

				} else {

					std::cerr << "Error: Wrong FASTA format: "
							<< sequenceFilepath_ << std::endl;
					exit( -1 );
				}
			}
		}

		if( !( header.empty() ) ){				// store the last sequence

			if( !( sequence.empty() ) ){

				size_t L = sequence.length();

				if ( L > maxL ){
					maxL = L;
				}
				if ( L < minL ){
					minL = L;
				}

				// translate sequence into encoding
				uint8_t* encoding = ( uint8_t* )calloc( L, sizeof( uint8_t ) );

				for( size_t i = 0; i < L; i++ ){

					encoding[i] = Alphabet::getCode( sequence[i] );

					if( encoding[i] == 0 ){

//						std::cerr << "Warning: The FASTA file contains an undefined base: " <<
//								sequence[i] << " at sequence " << header << std::endl;

						continue; // exclude undefined base from base counts
					}
					baseCounts[encoding[i]-1]++; // count base
				}
				sequences_.push_back( new Sequence( encoding, L, header, Y_, singleStrand ) );

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
		exit( -1 );
	}

	maxL_ = maxL;
	minL_ = minL;

	 // calculate the sum of bases
	size_t sumCounts = 0;
	for( size_t i = 0; i < Y_[1]; i++ ){
		sumCounts += baseCounts[i];
	}

	// calculate base frequencies
	for( size_t i = 0; i < Y_[1]; i++ ){
		baseFrequencies_[i] = static_cast<float>( baseCounts[i] ) /
		                      static_cast<float>( sumCounts );
	}

	return 0;
}

int SequenceSet::readIntensities(){

	std::cerr << "Error: SequenceSet::readIntensities() "
			"is not implemented so far." << std::endl;
	exit( -1 );
}
