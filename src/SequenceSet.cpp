/*
 * SequenceSet.cpp
 *
 *  Created on: Apr 7, 2016
 *      Author: wanwan
 */

#include <iostream>
#include <string>			/* string, find_last_of, memcpy */

#include <stdint.h>         /* uint8_t */

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

void createRevComp( SequenceSet* sequence_set ){
    // create reverse complement sequences for the sequence set
    for( uint8_t i = 1; i < sequence_set->N_; i++ ){
        unsigned int l = sequence_set->sequences_->L_;

        // reallocate the memory for storing reverse complement plus one after original sequence
        sequence_set->sequences_->sequence_ = ( char* )realloc( sequence_set->sequences_->sequence_, l*2 + 2);

        // set a 'N' character in the middle
        sequence_set->sequences_->sequence_[l+1] = 'N';
        for( int m = l+2, n = l; m <= 2*l + 1; m++, n-- ){
            switch( sequence_set->sequences_->sequence_[n] ){
            case 'A':  sequence_set->sequences_->sequence_[m] = 'T'; break;
            case 'C':  sequence_set->sequences_->sequence_[m] = 'G'; break;
            case 'G':  sequence_set->sequences_->sequence_[m] = 'C'; break;
            case 'T':  sequence_set->sequences_->sequence_[m] = 'A'; break;
            }
        }
    }
    // update sequenceSet information
    sequence_set->maxL_ =  sequence_set->maxL_*2 + 1;
    sequence_set->minL_ =  sequence_set->minL_*2 + 1;
}

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
	unsigned int sum = 0;										// sum count of all the mono-nucleotides
	unsigned int countA = 0, countT = 0, countC = 0, countG = 0;

	std::ifstream sequence_file( sequenceFilepath );			// read file from path
	std::string single_line;									// read each line from file
	std::string each_sequence;									// get each sequence
	std::vector<std::string> sequences;                         // store all the sequences in a vector array for sequence shuffling

	getline( sequence_file, each_sequence );
	if( each_sequence[0] != '>' ){								// read the first line, exit when it does not start with '>'
		fprintf( stderr, "file is not in FASTA format: "
		        "Header does not start with \">\" " );
		exit( -1 );
	} else {
		getline( sequence_file, each_sequence, '>' );			// read the first sequence on the second line
		sequences_->L_ = each_sequence.length();
		minL_ = maxL_ = sequences_->L_;						// initialize minL_ and maxL_
	}

	while( getline( sequence_file, single_line ).good() ){
		if( single_line[0] == '>' ){
			N_++;												// count the sequence number
			sequences_->header_() = single_line.substr( 1 );    // Take the header after the '>' sign.
			if( each_sequence != NULL ){						// compute from the second sequence on
			    sequences_->L_ = each_sequence.length();
				sum += sequences_->L_;                         // count the total number of all nucleotides
				if( sequences_->L_ > maxL_ ){
					maxL_ = sequences_->L_;
				}
				if( sequences_->L_ < minL_ ){
					minL_ = sequences_->L_;
				}
				// Either ...
				sequences.push_back( each_sequence );
				// Or ...
				sequences_->sequence_() = each_sequence;

				// count the occurrence of mono-nucleotides
				countA += std::count(each_sequence.begin(), each_sequence.end(), 'A');
				countC += std::count(each_sequence.begin(), each_sequence.end(), 'C');
				countG += std::count(each_sequence.begin(), each_sequence.end(), 'G');
				countT += std::count(each_sequence.begin(), each_sequence.end(), 'T');

				each_sequence.clear();
			}

		} else {
			each_sequence += single_line;
		}
	}

	/*
	 * Calculate the frequencies of all mono-nucleotides
	 * Then declare baseFrequencies_
	 */
	baseFrequencies_ = new float[4];
	baseFrequencies_ = { countA / (float)sum, countC / (float)sum, countG / (float)sum, countT / (float)sum };

	return 0;
}




