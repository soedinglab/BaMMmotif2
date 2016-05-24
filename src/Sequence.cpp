/*
 * Sequence.cpp
 *
 *  Created on: Apr 18, 2016
 *      Author: administrator
 */

#include <cstring>			// memcpy
#include <array>

#include "Sequence.h"
#include "Global.h"
#include "Alphabet.h"

Sequence::Sequence( uint8_t* sequence, int L, std::string header ){

	if( ! Global::revcomp ){
		L_ = L;
		sequence_ = ( uint8_t* )malloc( L_ );
		std::memcpy( sequence_, sequence, L_ );
	} else {
		L_ = 2 * L;
		createRevComp( sequence, L );
	}
	header_.assign( header );
}

Sequence::~Sequence(){
//	std::cout << "This is a destructor for Sequence class. " << std::endl;
}

uint8_t* Sequence::createRevComp( uint8_t* sequence, int L ){
	sequence_ = ( uint8_t* )malloc( 2* L );
	for( int i = 0; i < L; i++ ){
		sequence_[i] = sequence[i];											// the original sequence
		sequence_[2*L-i-1] = Alphabet::getComplementCode( sequence[i] );	// the reverse complementary
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
	int y = 0;
	for( int n = 0; n <= k; n++ )
		y += sequence_[i-k+n] * pow( Alphabet::getSize(), n );
	return y;
}

