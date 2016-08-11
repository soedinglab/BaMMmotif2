#include "Sequence.h"

Sequence::Sequence( uint8_t* sequence, int L, std::string header, std::vector<int> Y, bool revcomp ){

	if( revcomp ){
		L_ = 2 * L + 1;
		sequence_ = ( uint8_t* )calloc( L_, sizeof( uint8_t) );
		appendRevComp( sequence, L );
	} else{
		L_ = L;
		sequence_ = ( uint8_t* )calloc( L_, sizeof( uint8_t) );
		std::memcpy( sequence_, sequence, L_ );
	}
	header_ = header;
	Y_ = Y;
}

Sequence::~Sequence(){
	if( sequence_ != NULL ){
		free( sequence_ );
	}
}

uint8_t* Sequence::getSequence(){
	return sequence_;
}

int Sequence::getL(){
	return L_;
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
	 *  extract (k+1)mer y from positions (i-k,...,i) of the sequence
	 *  e.g.		| monomer |                      dimer                      |      trimer     ... | ...
	 *  (k+1)mer:	| A C G T | AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT | AAA AAC AAG AAT ... | ...
	 *  y:			| 0 1 2 3 |	 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 |  0   1   2   3  ... | ...
	 */

	int y = 0;
	for( int j = k; j >= 0; j-- ){
		if( sequence_[i-j] > 0 ){
			y += ( sequence_[i-j] -1 ) * Y_[j];
		} else{
			y = -1; // for non-defined alphabet letters
			break;
		}
	}
	return y;
}

void Sequence::print(){

	std::cout << ">" << header_ << std::endl;
	for( int i = 0; i < L_; i++ ){
		std::cout << Alphabet::getBase( sequence_[i] );
	}
	std::cout << std::endl;
}

uint8_t* Sequence::appendRevComp( uint8_t* sequence, int L ){

	for( int i = 0; i < L; i++ ){
		sequence_[i] = sequence[i];										// the sequence
		sequence_[2*L-i] = Alphabet::getComplementCode( sequence[i] );	// and its reverse complement
	}
	return sequence_;
}
