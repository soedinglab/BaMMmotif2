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

Sequence::Sequence( uint8_t* sequence, unsigned int LS, std::string header ){

	if( ! Global::revcomp ){
		std::memcpy( sequence_, sequence, LS );
	} else {
		createRevComp( sequence );
	}

	L_ = LS;
	header_.assign( header );

}

Sequence::~Sequence(){
	std::cout << " This is a distructor to be fulfilled. " << std::endl;
}

uint8_t* Sequence::createRevComp( uint8_t* seq ){
	for( unsigned int i = 0; i < L_; i++ ){								// deep copy
		sequence_[i] = seq[i];											// the original sequence
		sequence_[2*L_-i -1] = Alphabet::getComplementCode( seq[i] );	// the reverse complementary
	}
	return sequence_;
}

unsigned int Sequence::getL(){
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

