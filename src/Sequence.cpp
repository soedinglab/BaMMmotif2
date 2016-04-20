/*
 * Sequence.cpp
 *
 *  Created on: Apr 18, 2016
 *      Author: administrator
 */

#include <cstring>			// memcpy

#include <Sequence.h>

Sequence::Sequence( uint8_t* sequence, int L, char* header, float intensity = 0.0f, float weight = 1.0f ){
	if( ! Global::revcomp ){
		std::memcpy( sequence_, sequence, L ) ;
	} else {
		createRevComp( sequence );
	}

	L_ = L;
	header_ = header;
	setIntensity( intensity );
	setWeight( weight );
}

uint8_t* Sequence::createRevComp( uint8_t* seq ){
	for( int i = 0; i < L_; i++ ){										// deep copy
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

char* Sequence::getHeader(){
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
