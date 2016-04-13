/*
 * SequenceSet.cpp
 *
 *  Created on: Apr 7, 2016
 *      Author: wanwan
 */

#include <iostream>
#include <string>			/* string, find_last_of, memcpy */

#include "SequenceSet.h"

char* extractFilenameFromPath( char* filePath ){
    char* fileName;
    std::string path = filePath;
    std::size_t lastSlashIndex = path.find_last_of( "/\\" );
    if( std::string::npos != lastSlashIndex ){
        path.erase( 0, lastSlashIndex + 1 );
    }
    std::size_t periodIndex = path.rfind( '.' );
    if( std::string::npos != periodIndex ){
        path.erase( periodIndex );
    }
    fileName = path;            // after removing directory and extension
    return fileName;
};

int SequenceSet::readFASTA( char* sequenceFilepath ){

	/*
	 * while reading in the sequences, do:
	 * 1. extract the header for each sequence
	 * 2. calculate the minimal and maximal length in the file
	 * 3. count the number of sequences
	 * 4. calculate mono-nucleotide frequencies
	 */
    filename_ = extractFilenameFromPath( sequenceFilepath );
	N_ = 0;
	unsigned int sequence_length;
	unsigned int sIndex;
	unsigned int sum = 0;										// sum count of all the mono-nucleotides
	unsigned int countA = 0, countT = 0, countC = 0, countG = 0;
	float freqA, freqT, freqC, freqG;

	std::ifstream sFile( sequenceFilepath );					// read file from path
	std::string sLine;											// read each line from file
	std::string eachSequence;									// get each sequence
	std::vector<std::string> sequences;

	getline( sFile, eachSequence );
	if( eachSequence[0] != '>' ){								// read the first line, exit when it does not start with '>'
		fprintf( stderr, "file is not in FASTA format: "
		        "Header does not start with \">\" " );
		exit( -1 );
	} else {
		getline( sFile, eachSequence, '>' );					// read the first sequence on the second line
		sequence_length = eachSequence.length();
		minL_ = maxL_ = sequence_length;						// initialize minL_ and maxL_
	}

//	sIndex = sequence_length;
//	char* sequences;											// a dynamic array to store sequences
//	const char* sequence = ( sIndex );							// a const array which has a fixed size of the first sequence

	while( getline( sFile, sLine ).good() ){
		if( sLine[0] == '>' ){
			N_++;												// count the sequence number
			sequences_->header_() = sLine.substr( 1 );          // Take the header after the '>' sign.
			if( eachSequence != NULL ){							// compute from the second sequence on
				sequence_length = eachSequence.length();
				sum += sequence_length;                         // count the total number of all nucleotides
				if( sequence_length > maxL_ ){
					maxL_ = sequence_length;
				}
				if( sequence_length < minL_ ){
					minL_ = sequence_length;
				}
				sequences.push_back( eachSequence );            // Either ...
				sequences_->sequence_() = eachSequence;         // Or ...

				// count the occurrence of mono-nucleotides
				countA += std::count(eachSequence.begin(), eachSequence.end(), 'A');
				countC += std::count(eachSequence.begin(), eachSequence.end(), 'C');
				countG += std::count(eachSequence.begin(), eachSequence.end(), 'G');
				countT += std::count(eachSequence.begin(), eachSequence.end(), 'T');

				eachSequence.clear();
			}

		} else {
			eachSequence += sLine;
		}
	}

	/*
	 * Calculate the frequencies of mono-nucleotides
	 */
	baseFrequencies_[0] = freqA = (float)countA / (float)sum;
	baseFrequencies_[1] = freqC = (float)countC / (float)sum;
	baseFrequencies_[2] = freqG = (float)countG / (float)sum;
	baseFrequencies_[3] = freqT = (float)countT / (float)sum;

	return 0;
}



