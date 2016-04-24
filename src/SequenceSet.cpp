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


SequenceSet::SequenceSet( char* sequenceFilepath, char* intensityFilepath ){

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

std::vector<Sequence*> SequenceSet::getSequences(){
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
	unsigned int LS;								// length of the sequence
	unsigned int LH;								// length of the header
	unsigned int sizeSeq = 1000;					// initialize the size of sequence array as 1000
	unsigned int* baseCounts = ( unsigned int* )calloc( Alphabet::getSize() + 1, sizeof( unsigned int ) );
	unsigned int baseSumCount = 0;					// total number of nucleotides in the FASTA file
	unsigned int maxL = 0, minL = 1000;				// initialize the min. and max. length

	std::string header;
	uint8_t* sequence = ( uint8_t* )calloc( sizeSeq, sizeof( uint8_t ) );
	FILE* fp;
	Sequence* seq;

	if( ( fp = fopen( sequenceFilepath_, "r" ) ) == NULL ){
		fprintf( stderr, "Cannot open positive sequence file: %s\n", sequenceFilepath_ );
		exit( -1 );
	}

	int base; // in ASCII code
	base = fgetc( fp );

	fprintf( stderr, "Start reading in the fasta file! \n" );

	// skipping leading blank lines
	while( ( base == '\n' || base == '\r' ) && base != EOF ){
		base = fgetc( fp );
	}

	while( base != EOF ){
		// reset the parameters to 0 or empty
		LS = 0;
		LH = 0;
		header.clear();
		memset( sequence, 0, sizeSeq );

		if( base != '>' ){
			fprintf( stderr, "Error: Not a FASTA format file! \n"
					"Please check the content in the file. \n" );
			exit( -1 );
		}

		// if everything is fine, now base = '>'

		// count the number of header
		N++;

		base = fgetc( fp ); // read in the first char after '>' in the header
		fprintf( stderr, "Start reading in the header! \n" );

		while( ( base != '\n' && base != '\r' ) && base != EOF ) {
			LH++;
			// append each letter to the header
			header += ( char )base;
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
			// count the length of the sequence
			LS++;

			// count the occurrence of mono-nucleotides and sum them up
			baseCounts[ Alphabet::getCode( ( char )base ) ]++;				// convert char to int

			// before memcpy, check if the sequence length exceeds the size of array
			if( LS > sizeSeq ){
				sizeSeq *= 2;
				sequence = ( uint8_t* )realloc( sequence, sizeSeq );
			}

			// memcpy sequence into sequence array
			std::memcpy( &sequence[LS-1], &base, 1 );

			// read in the next nucleobase
			base = fgetc( fp );
		}

		fprintf( stderr, "Finish reading in the sequence! \n" );
		fprintf( stderr, "The header is %s \n", header.c_str() );
		fprintf( stderr, "The sequence is %s \n", sequence );
		fprintf( stderr, "The length of sequence is %d \n", LS );
		for(int i = 1; i < 5; i++ ){
			fprintf( stderr, " The base count is %d \n", baseCounts[i] );
		}
		// create the sequence object to sequences_
		seq = new Sequence( sequence, LS, header );
		fprintf( stderr, "A new sequence object is created! \n" );

		sequences_.push_back( seq );

		fprintf( stderr, "A new sequence object is push to the sequences vector! \n" );

		// check if sequence is empty
		if( LS == 0 ){
			fprintf( stderr, "Empty sequence with the entry %d !\n ", N );
		}

		// check the length of each sequence if it is min or max
		if ( LS > maxL )
			maxL = LS;
		if ( LS < minL )
			minL = LS;
		fprintf( stderr, "header of sequence is %s \n", sequences_[N-1]->getHeader().c_str());
//		fprintf( stderr, "the sequence is %s \n", sequences_[N-1]->getSequence() );
//		fprintf( stderr, "the length is %d \n", sequences_[N-1]->getL() );
	}

	for( unsigned int i = 0; i < N; i++ ){
		fprintf( stderr, "header of sequence is %s \n", sequences_[i]->getHeader().c_str());
		fprintf( stderr, "the sequence is %s \n", sequences_[i]->getSequence() );
		fprintf( stderr, "the length is %d \n", sequences_[i]->getL() );
	}
	fprintf( stderr, "I am here! \n" );
	exit( -1 );

	// Calculate the frequencies of all mono-nucleotides
	for( unsigned int i = 0; i < Alphabet::getSize(); i++ ){
		baseSumCount += baseCounts[i+1];
	}
	for( unsigned int i = 0; i < Alphabet::getSize(); i++ ){
		baseFrequencies_[i] = ( float )baseCounts[i+1] / ( float )baseSumCount;
	}

	// assign the maxL_, minL_ and N_
	maxL_ = maxL;
	minL_ = minL;
	N_ = N;

	delete baseCounts;
	delete sequence;
	delete fp;
	delete seq;

	fprintf( stderr, "max length is %d \n", maxL_ );
	fprintf( stderr, "min length is %d \n", minL_ );
	fprintf( stderr, "number of sequences is %d \n", N_ );
	fprintf( stderr, "I am here! \n" );
	exit( -1 );

	return 0;
}

int SequenceSet::readIntensities(){
	return -1;
}

SequenceSet::~SequenceSet(){
	std::cout << " This is a distructor to be fulfilled. " << std::endl;
	delete sequenceFilepath_;
	if( intensityFilepath_ ) delete intensityFilepath_;
//	delete sequences_;
//	delete baseFrequencies_;
}
