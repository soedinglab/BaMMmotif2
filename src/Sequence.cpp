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
		L_ = 2 * L + 1;
		sequence_ = ( uint8_t* )malloc( L_ );
		createRevComp( sequence, L );
	}
	header_.assign( header );
}

Sequence::~Sequence(){
//	std::cout << "This is a destructor for Sequence class. " << std::endl;
//	if( sequence_ ) free( sequence_ );
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
	int y = 0;
	if( sequence_[i] > 0 )
		for( int n = 0; n <= k; n++ )
			y += ( sequence_[i-k+n] -1 ) * Global::ipow( Alphabet::getSize(), n );
	else
		y = -1;
	return y;
}

int	Sequence::extractKmerbg( int i, int k ){
	int y = -1;
	if( sequence_[i] != 0 )
		for( int n = 0; n <= k; n++ )
			y += sequence_[i-k+n] * Global::ipow( Alphabet::getSize(), n );
	return y;
}
