#include "Sequence.h"
#include "../refinement/utils.h"

Sequence::Sequence( uint8_t* sequence,
					size_t L,
					std::string header,
					std::vector<size_t> Y,
					bool singleStrand ){

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

	for( size_t i = 0; i < 12; i++ ){
		Y_.push_back( ipow( Alphabet::getSize(), i ) );
	}

	/**
	 *  extract (k+1)-mer y from positions (i-k,...,i) of the sequence
	 *  e.g.		| monomer |                      dimer                      |      trimer     ... | ...
	 *  (k+1)mer:	| A C G T | AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT | AAA AAC AAG AAT ... | ...
	 *  y:			| 0 1 2 3 |	 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 |  0   1   2   3  ... | ...
	 *
	 ** note:
	 *  randomize the unknown letter N to A, C, G or T.
	 */

	kmer_ = ( size_t* )calloc( L_, sizeof( size_t ) );
	for( size_t i = 0; i < L_; i++ ){
		for( size_t k = i < 10 ? i+1 : 11; k > 0; k-- ){
			kmer_[i] += ( ( sequence_[i-k+1] == 0 ) ? ( size_t )rand() % Y_[1] :
						( sequence_[i-k+1] - 1 ) ) * Y_[k-1];
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
//		std::cout << kmer_[i] << '\t';
	}
	std::cout << std::endl;
}

void Sequence::appendRevComp( uint8_t* sequence, size_t L ){

	for( size_t i = 0; i < L; i++ ){
		// the forward sequence
		sequence_[i] = sequence[i];
		// the reverse complement
		sequence_[2*L-i] = Alphabet::getComplementCode( sequence[i] );
	}
}
