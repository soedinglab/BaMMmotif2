/*
 * Sequence.cpp
 *
 *  Created on: Apr 18, 2016
 *      Author: administrator
 */

#include <cstring>			// memcpy
#include <array>

#include "Sequence.h"

#ifndef VIRUS_X_
#include "../bamm/Global.h"
#else
#include "../virusX/Global.h"
#endif

#include "Alphabet.h"

Sequence::Sequence( uint8_t* sequence, int L, std::string header ){

	if( !Global::revcomp ){
		L_ = L;
//		sequence_ = sequence;
		sequence_ = ( uint8_t* )malloc( L_ );
		std::memcpy( sequence_, sequence, L_ );
	} else {
		L_ = 2 * L + 1;
		sequence_ = ( uint8_t* )malloc( L_ );
		createRevComp( sequence, L );
	}
	header_.assign( header );
}

//Sequence::Sequence( const Sequence& other ){
//	L_ = other.L_;
//	sequence_ = other.sequence_;
//	header_ = other.header_;
//}

Sequence::~Sequence(){
//	if( sequence_ ) free( sequence_ );
//	std::cout << "Destructor for Sequence class works fine. \n";			// it will be print thousands of times
}

uint8_t* Sequence::createRevComp( uint8_t* sequence, int L ){

	for( int i = 0; i < L; i++ ){
		sequence_[i] = sequence[i];											// the original sequence
		sequence_[i+1] = 0;													// add a '0' at the junction site
		sequence_[2*L-i] = Alphabet::getComplementCode( sequence[i] );		// the reverse complementary
	}
	return sequence_;
}

int Sequence::getL(){
	return L_;
}

uint8_t* Sequence::getSequence(){
	return sequence_;
}

std::string Sequence::getHeader(){
	return header_;
}

float Sequence::getIntensity(){
	return intensity_;
}

float Sequence::getWeight(){
	return weight_;
}

void Sequence::setIntensity( float intensity ){
	intensity_ = intensity;
}

void Sequence::setWeight( float weight ){
	weight_ = weight;
}

int	Sequence::extractKmer( int i, int k ){

	/**
	 *  extract (k+1)-mer y from positions (i-k,...,i) of the sequence
	 *  e.g.			|  mono-mer	 |								dimer							  |		trimer		...	| ...
	 *  (k+1)-mer:		A	C	G	T	AA	AC	AG	AT	CA	CC	CG	CT	GA	GC	GG	GT	TA	TC	TG	TT	AAA	AAC	AAG	AAT	...
	 *  extracted y:	0	1	2	3|	0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15|	0	1	2	3	...
	 */

	int n, y = 0;
	for( n = k; n >= 0; n-- ){
		if( sequence_[i-n] > 0 )
			y += ( sequence_[i-n] -1 ) * Global::powA[n];
		else {
			y = -1;										// for unknown alphabets, set y to -1
			break;
		}
	}
	return y;
}
