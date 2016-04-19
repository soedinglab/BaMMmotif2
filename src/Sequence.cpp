/*
 * Sequence.cpp
 *
 *  Created on: Apr 18, 2016
 *      Author: administrator
 */

#include <Sequence.h>

Sequence::Sequence( uint8_t* sequence, int L, char* header, float intensity, float weight ){
	sequence_ = sequence ;
	L_ = L;
	header_ = header;
	setIntensity( intensity );
	setWeight( weight );
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
