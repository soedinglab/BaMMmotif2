/*
 * SequenceSet.cpp
 *
 *  Created on: Apr 7, 2016
 *      Author: wanwan
 */

#include <iostream>			// std::
#include <fstream>			// std::fstream
#include <cstdio>           // EOF
#include <string>			// string, find_last_of,
#include <cstring>			// memcpy
#include <algorithm>		// count

#include "SequenceSet.h"
#include "Alphabet.h"


SequenceSet::SequenceSet( char* sequenceFilepath, char* intensityFilepath ){

	sequenceFilepath_ = sequenceFilepath;
	intensityFilepath_ = intensityFilepath;
	baseFrequencies_ = new float[Alphabet::getSize()];
//	sequences_ = ( Sequence* )calloc( 5000, sizeof( Sequence ) );

	readFASTA();

	// if intensity file is given, read in intensity file
	if ( intensityFilepath != NULL )
		readIntensities();
}

SequenceSet::~SequenceSet(){
	if( baseFrequencies_ ) delete []baseFrequencies_;
	std::cout << "Destructor for SequenceSet class works fine. \n";
}

char* SequenceSet::getSequenceFilepath(){
	return sequenceFilepath_;
}

char* SequenceSet::getIntensityFilepath(){
	return intensityFilepath_;
}

std::vector<Sequence> SequenceSet::getSequences(){
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

int SequenceSet::readFASTA(){

	/**
	 * while reading in the sequences, do:
	 * 1. extract the header for each sequence
	 * 2. calculate the minimal and maximal length in the file
	 * 3. count the number of sequences
	 * 4. calculate mono-nucleotide frequencies
	 */
/*
	unsigned int N = 0; 							// sequence counter
	unsigned int ls;								// length of the sequence
	unsigned int lh;								// length of the header
	unsigned int sizeSeq = 1000;					// initialize the size of sequence array as 1000
//	unsigned int* baseCounts = ( unsigned int* )calloc( Alphabet::getSize() , sizeof( unsigned int ) );
	unsigned int* baseCounts = new unsigned int[Alphabet::getSize()];
	unsigned int baseSumCount = 0;					// total number of nucleotides in the FASTA file
	unsigned int maxL = 0, minL = 1000;				// initialize the min. and max. length

	FILE* fp;
	std::string header;

//	int* sequence = ( int* )calloc( sizeSeq, sizeof( int ) );
	int* sequence = new int[sizeSeq];

	if( ( fp = fopen( sequenceFilepath_, "r" ) ) == NULL ){
		fprintf( stderr, "Cannot open positive sequence file: %s\n", sequenceFilepath_ );
		exit( -1 );
	}

	int base; // in ASCII code
	base = fgetc( fp );

	fprintf( stderr, "Start reading in the FASTA file! \n" );

	// skipping leading blank lines
	while( ( base == '\n' || base == '\r' ) && base != EOF ){
		base = fgetc( fp );
	}

	while( base != EOF ){

		// reset parameters to 0
		ls = 0;
		lh = 0;
		header = "";
		memset( &sequence, 0, ls );

		if( base != '>' ){
			fprintf( stderr, "Error: Not a FASTA format file! \n"
					"Please check the content in the file. \n" );
			exit( -1 );
		}

		N++;										// count the number of header

		base = fgetc( fp ); 						// read in the first char after '>' in the header
		fprintf( stderr, "Start reading in each header... \n" );

		while( ( base != '\n' && base != '\r' ) && base != EOF ) {
			lh++;
			header += ( char )base;					// append each letter to the header
			base = fgetc( fp );
		} 											// end of header or start of sequence line

		while( base == '\n' || base == '\r' ){		// skipping blank lines in between
			base = fgetc( fp );
		} 											// start of sequence line

		fprintf( stderr, "Start reading in each sequence... \n" );
		while( base != '>' && base != EOF ){

			if( base == '\n' || base == '\r' ){		// skip all the '\n' and '\r'
				base = fgetc( fp );
				continue;
			}

			ls++;

			baseCounts[ Alphabet::getCode( ( char )base )-1 ]++;// count the occurrence of mono-nucleotides

			if( ls > sizeSeq ){									// before memcpy, check if the sequence length exceeds the size of array
				sizeSeq *= 2;
				sequence = ( int* )realloc( sequence, sizeSeq );
			}

//			std::memcpy( &sequence[ls-1], &base, 1 );			// memcpy base into sequence array
			strncat( sequence, base, 1 );
		}

		Sequence* seq = new Sequence( sequence, ls, header );
		fprintf( stderr, "A new sequence object is created! \n" );

		sequences_.push_back( seq );
		fprintf( stderr, "A new sequence object is push to the sequences vector! \n" );

		if( ls == 0 ){											// check if sequence is empty
			fprintf( stderr, "Empty sequence with the entry %d !\n ", N );
		}

		if ( ls > maxL )										// check the length of each sequence if it is min or max
			maxL = ls;
		if ( ls < minL )
			minL = ls;

		fprintf( stderr, "header of sequence is %s \n", sequences_[N-1]->getHeader().c_str());
		fprintf( stderr, "the length is %d \n", sequences_[N-1]->getL() );
		for( int i = 1; i < 5; i++ ){
			fprintf( stderr, " The base count is %d \n", baseCounts[i] );
		}

		delete seq;
		base = fgetc( fp );

	}

	fprintf( stderr, "Finish reading in the FASTA file. \n" );


	for( unsigned int i = 0; i < Alphabet::getSize(); i++ ){	// Calculate the sum of all mono-nucleotides
		baseSumCount += baseCounts[i+1];
	}
	for( unsigned int i = 0; i < Alphabet::getSize(); i++ ){	// Calculate the frequencies of all mono-nucleotides
		baseFrequencies_[i] = ( float )baseCounts[i+1] / ( float )baseSumCount;
	}

	// assign the maxL_, minL_ and N_
	maxL_ = maxL;
	minL_ = minL;
	N_ = N;

	free( sequence );
	delete [] baseCounts;

	fprintf( stderr, "max length is %d \n", maxL_ );
	fprintf( stderr, "min length is %d \n", minL_ );
	fprintf( stderr, "number of sequences is %d \n", N_ );
*/

	std::ifstream seqfile( sequenceFilepath_ ); 				// open sequence file
	std::string line;
	std::string header;
	std::string sequence;

	int maxL = 0, minL = 1000;									// initialize the min. and max. length
	int N = 0;													// count sequences
	int i;
	unsigned int* baseCounts = new unsigned int[Alphabet::getSize()];
	for( i = 0; i < Alphabet::getSize(); i++ )
		baseCounts[i] = 0;

	if( seqfile ){
		while( getline( seqfile, line ).good() ){
			if( line.empty() || line[0] == '>' ){
				if( !header.empty() ){
					int L = ( int )sequence.length();
					uint8_t* encodeSeq = new uint8_t[L];

					for( i = 0; i < L; i++ ){
						encodeSeq[i] = Alphabet::getCode( sequence[i] );
						if( encodeSeq[i] == 0 ){
							fprintf( stderr, "Warning: The FASTA file contains other alphabet(s) "
									"than %s on line %d.\n", Alphabet::getAlphabet(), N );
							continue;
						}
						baseCounts[encodeSeq[i]-1]++;			// count the occurrence of mono-nucleotides
					}

					Sequence seq( encodeSeq, L, header );
					sequences_.push_back( seq );

					if ( L > maxL )								// check the length of each sequence if it is min or max
						maxL = L;
					if ( L < minL )
						minL = L;

					N++;

					sequence.clear();
					header.clear();
					delete []encodeSeq;
				}
				if( !line.empty()  ){
					if( line.length() == 1 )					// for empty entry
						header = "NonName";
					else
						header = line.substr(1);				// fetch the header
				}
			} else if( !header.empty() ){
				if( line.find(' ') != std::string::npos ){		// Invalid sequence--no spaces allowed
					header.clear();
					sequence.clear();
				} else {
					sequence += line;
				}
			}
		}
		if( !header.empty() ){									// for the last entry if no new line follows by
			int L = ( int )sequence.length();
			uint8_t* encodeSeq = new uint8_t[L];
			for( int i = 0; i < L; i++ ){
				encodeSeq[i] = Alphabet::getCode( sequence[i] );
				if( encodeSeq[i] == 0 ){
					fprintf( stderr, "Warning: The FASTA file contains other alphabet(s) "
							"than %s on line %d.\n", Alphabet::getAlphabet(), N );
					continue;
				}
				baseCounts[encodeSeq[i]-1]++;					// count the occurrence of mono-nucleotides
			}
			Sequence seq( encodeSeq, L, header );
			sequences_.push_back( seq );

			if ( L > maxL )										// check the length of each sequence if it is min or max
				maxL = L;
			if ( L < minL )
				minL = L;

			N++;

			sequence.clear();
			header.clear();
			delete []encodeSeq;
		}
	} else {
		fprintf( stderr, "Error: Cannot open positive sequence file: %s\n", sequenceFilepath_ );
		exit( -1 );
	}
	seqfile.close();											// close sequence file

	if( minL == 0 ){											// exclude format with empty sequence under '>'
		fprintf( stderr, "Error: Wrong FASTA format file!\n" );
		exit( -1 );
	}

	maxL_ = maxL;
	minL_ = minL;
	N_ = N;

	unsigned int baseSumCount = 0;
	for( i = 0; i < Alphabet::getSize(); i++ )				// Calculate the sum of all mono-nucleotides
		baseSumCount += baseCounts[i];

	for( i = 0; i < Alphabet::getSize(); i++ )				// Calculate the frequencies of all mono-nucleotides
		baseFrequencies_[i] = ( float )baseCounts[i] / ( float )baseSumCount;

	delete []baseCounts;

	return 0;
}

int SequenceSet::readIntensities(){
	return -1;
}
