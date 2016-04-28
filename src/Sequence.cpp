/*
 * Sequence.cpp
 *
 *  Created on: Apr 18, 2016
 *      Author: administrator
 */

#include <cstring>			// memcpy

#include "Sequence.h"
#include "Global.h"
#include "Alphabet.h"

Sequence::Sequence( uint8_t* sequence, int L, std::string header ){

	if( ! Global::revcomp ){
		sequence_ = ( uint8_t* )malloc( L );
		std::memcpy( sequence_, sequence, L );
	} else {
		createRevComp( sequence );
	}

	L_ = L;
	header_.assign( header );
}

Sequence::~Sequence(){
//	std::cout << "This is a destructor for Sequence class. " << std::endl;
	// No need to free sequence_, => temporary objects are automatically destroyed when loop out
}

uint8_t* Sequence::createRevComp( uint8_t* sequence ){
	for( int i = 0; i < L_; i++ ){										// deep copy
		sequence_[i] = sequence[i];											// the original sequence
		sequence_[2*L_-i-1] = Alphabet::getComplementCode( sequence[i] );	// the reverse complementary
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

