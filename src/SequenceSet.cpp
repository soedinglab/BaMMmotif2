/*
 * SequenceSet.cpp
 *
 *  Created on: Apr 7, 2016
 *      Author: wanwan
 */

#include <iostream>			// std::
#include <fstream>			// std::fstream
#include <string>			// string, find_last_of,
#include <cstring>			// memcpy
#include <algorithm>		// std::count

#include "SequenceSet.h"
#include "Alphabet.h"

SequenceSet::SequenceSet( char* sequenceFilepath, char* intensityFilepath = NULL ){

	sequenceFilepath_ = sequenceFilepath;
	intensityFilepath_ = intensityFilepath;

	// read in FASTA file
	readFASTA();

	// calculate mono-nucleotide frequencies
	getBaseFrequencies();

	// if intensity file is given, read in intensity file
	if ( intensityFilepath != NULL )
		readIntensities();
}

char* SequenceSet::getSequenceFilepath(){
	return sequenceFilepath_;
}

char* SequenceSet::getIntensityFilepath(){
	return intensityFilepath_;
}

Sequence* SequenceSet::getSequences(){
	return sequences_;
}

unsigned int SequenceSet::getN(){
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

int SequenceSet::readFASTA(){

	/**
	 * while reading in the sequences, do:
	 * 1. extract the header for each sequence
	 * 2. calculate the minimal and maximal length in the file
	 * 3. count the number of sequences
	 * 4. calculate mono-nucleotide frequencies
	 */

	unsigned int N = 0; 							// sequence counter
	unsigned int LS = 0;							// length of the sequence
	unsigned int LH = 0;							// length of the header
	unsigned int sizeHeader = 100;					// initialize the size of header array as 100
	unsigned int sizeSeq = 1000;					// initialize the size of sequence array as 1000
	unsigned int* baseCounts = ( unsigned int* )calloc( Alphabet::getSize() + 1, sizeof( unsigned int ) );
	unsigned int baseSumCount = 0;					// total number of nucleotides in the FASTA file
	unsigned int maxL = 0, minL = 1000;				// initialize the min. and max. length

	char* header = ( char* )calloc( sizeHeader, sizeof( char ) );
	int* sequence = ( int* )calloc( sizeSeq, sizeof( int ) );

	FILE* fp;
	if( ( fp = fopen( sequenceFilepath_, "r" ) ) == NULL ){
		fprintf( stderr, "Cannot open positive sequence file: %s\n", sequenceFilepath_ );
		exit( -1 );
	}

	int base; // in ASCII code

	base = fgetc( fp );

	// skipping leading blank lines
	while( ( base == '\n' || base == '\r' ) && base != EOF ){
		base = fgetc( fp );
	}

	while( base != EOF ){
		// reset the length of sequence to 0
		LS = 0;
		LH = 0;

		if( base != '>' ){
			fprintf( stderr, "Error: Not a FASTA format file!" );
			exit( -1 );
		}

		// if everything is fine, now base = '>'

		// count the number of header
		N++;

		base = fgetc( fp ); // read in the first char after '>' in the header
		while( ( base != '\n' && base != '\r' ) && base != EOF ) {

			LH++;

			// memcpy header into header array
			std::memcpy( &header[LH-1], &base, 1);

			// what if the header will be empty? => copy '\n' into the header? => yes!
			base = fgetc( fp );
		} // end of header or start of sequence line

		// skipping blank lines in between
		while( base == '\n' || base == '\r' ){
			base = fgetc( fp );
		} // start of sequence line

		while( base != '>' && base != EOF ){

			// skip all the '\n' and '\r'
			if( base == '\n' || base == '\r' ){
				base = fgetc( fp );
				continue;
			}

			// if everything is fine, now base is the first nucleobase

			// count the length of the sequence
			LS++;

			// count the occurrence of mono-nucleotides and sum them up
			baseCounts[ Alphabet::getCode( base ) ]++;

			// before memcpy, check if the sequence length exceeds the size of array
			if( LS > sizeSeq ){
				sizeSeq *= 2;
				sequence = ( int* )realloc( sequence, sizeSeq );
			}

			// memcpy sequence into sequence array
			std::memcpy( &sequence[LS-1], &base, 1 );

			// read in the next nucleobase
			base = fgetc( fp );
		}

		// create the sequence object to sequences_
		sequences_[N-1] = new Sequence( sequence, LS, header );

		// check if sequence is empty
		if( LS == 0 ){
			fprintf( stderr, "Empty sequence with the entry %d !\n ", N );
		}

		// check the length of each sequence if it is min or max
		if ( LS > maxL )
			maxL = LS;
		if ( LS < minL )
			minL = LS;

	}

	// Calculate the frequencies of all mono-nucleotides
	for( unsigned int i = 0; i < Alphabet::getSize(); i++ ){
		baseSumCount += baseCounts[i+1];
	}
	for( unsigned int i = 0; i < Alphabet::getSize(); i++ ){
		baseFrequencies_[i] = baseCounts[i+1] / baseSumCount;
	}

	// assign the maxL_, minL_ and N_
	maxL_ = maxL;
	minL_ = minL;
	N_ = N;

/*
	std::ifstream sequenceFile( sequenceFilepath_ );			// read file from path
	if( !sequenceFile.is_open() ){
		fprintf( stderr, "Failed to open the sequence file!");
		exit( -1 );
	}
	std::string singleLine;										// read line after line from file
	std::string eachSequence;									// temporarily store each sequence
	std::vector<Sequence> sequences;							// store all the sequences in a vector, convenient for sequence shuffling
	std::string header;
	unsigned int maxL, minL;

	std::getline( sequenceFile, eachSequence );
	if( eachSequence[0] != '>' ){								// read the first line, exit when it does not start with '>'
		fprintf( stderr, "file is not in FASTA format: "
		        "Header does not start with \">\" " );
		exit( -1 );
	} else {
		header = eachSequence.substr( 1 );
		N++;
		std::getline( sequenceFile, eachSequence, '>' );		// read the first sequence
		minL = maxL = eachSequence.length();					// initialize minL and maxL
	}

	while( std::getline( sequenceFile, singleLine ).good() ){
		if( singleLine[0] == '>' ){
			Sequence sequence = new Sequence(eachSequence, eachSequence.length(), header );
			sequences.push_back( sequence );
			N++;
			header = singleLine.substr( 1 );

			if( sequence.getL() > maxL ){
				maxL = sequence.getL();
			}
			if( sequence.getL() < minL ){
				minL = sequence.getL();
			}

			// count the occurrence of mono-nucleotides
			for( int i = 0; i < Alphabet::getSize(); i++ ){
				baseCounts[i+1] += std::count(eachSequence.begin(), eachSequence.end(), Alphabet::getAlphabet()[i] );
				baseSumCount += baseCounts[i+1];
			}
			eachSequence.clear();

		} else {
			eachSequence += singleLine;
		}
	}
	Sequence sequence = new Sequence(eachSequence, eachSequence.length(), header );
	sequences.push_back( sequence );

	N_ = N;
	sequences_ = sequences;
	maxL_ = maxL;
	minL_ = minL;
*/

	return 0;
}




