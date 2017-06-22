#include "Sequence.h"
#include "utils.h"

Sequence::Sequence( uint8_t* sequence, size_t L, std::string header, std::vector<size_t> Y, bool singleStrand ){

	if( !singleStrand ){
		L_ = 2 * L + 1;
		sequence_ = ( uint8_t* )calloc( L_, sizeof( uint8_t ) );
		appendRevComp( sequence, L );
	} else {
		L_ = L;
		sequence_ = ( uint8_t* )calloc( L_, sizeof( uint8_t ) );
		std::memcpy( sequence_, sequence, L_ );
	}
	header_ = header;

	for( size_t i = 0; i <= 11; i++ ){
		Y_.push_back( ipow( Alphabet::getSize(), i ) );
	}

	/**
	 *  extract (k+1)mer y from positions (i-k,...,i) of the sequence
	 *  e.g.		| monomer |                      dimer                      |      trimer     ... | ...
	 *  (k+1)mer:	| A C G T | AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT | AAA AAC AAG AAT ... | ...
	 *  y:			| 0 1 2 3 |	 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 |  0   1   2   3  ... | ...
	 */

	kmer_ = ( int* )calloc( L_, sizeof( int ) );
	for( size_t i = 0; i < L_; i++ ){
		for( size_t j = i < 8 ? ( i+1 ) : 9; j > 0; j-- ){
			kmer_[i] += ( static_cast<int>( sequence_[i-j+1] ) - 1 ) * static_cast<int>( Y_[j-1] );
		}
	}
}

Sequence::~Sequence(){
	if( sequence_ != NULL ){
		free( sequence_ );
	}

	free( kmer_ );
}

uint8_t* Sequence::getSequence(){
	return sequence_;
}

size_t Sequence::getL(){
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

void Sequence::print(){

	std::cout << ">" << header_ << std::endl;
	for( size_t i = 0; i < L_; i++ ){
		std::cout << Alphabet::getBase( sequence_[i] );
	}
	std::cout << std::endl;
}

void Sequence::appendRevComp( uint8_t* sequence, size_t L ){

	for( size_t i = 0; i < L; i++ ){
		sequence_[i] = sequence[i];										// the sequence
		sequence_[2*L-i] = Alphabet::getComplementCode( sequence[i] );	// and its reverse complement
	}
}
