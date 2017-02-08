#include "Sequence.h"
#include "utils.h"

Sequence::Sequence( uint8_t* sequence, int L, std::string header, std::vector<int> Y, bool revcomp ){

	if( revcomp ){
		L_ = 2 * L + 1;
		sequence_ = ( uint8_t* )calloc( L_, sizeof( uint8_t ) );
		appendRevComp( sequence, L );
	} else {
		L_ = L;
		sequence_ = ( uint8_t* )calloc( L_, sizeof( uint8_t ) );
		std::memcpy( sequence_, sequence, L_ );
	}
	header_ = header;
	Y_ = new int[10];
	for( int i = 0; i < 10; i++ ){
		Y_[i] = ipow( Alphabet::getSize(), i );
	}

}

Sequence::~Sequence(){
	if( sequence_ != NULL ){
		free( sequence_ );
	}
	delete[] Y_;
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

void Sequence::print(){

	std::cout << ">" << header_ << std::endl;
	for( int i = 0; i < L_; i++ ){
		std::cout << Alphabet::getBase( sequence_[i] );
	}
	std::cout << std::endl;
}

void Sequence::appendRevComp( uint8_t* sequence, int L ){

	for( int i = 0; i < L; i++ ){
		sequence_[i] = sequence[i];										// the sequence
		sequence_[2*L-i] = Alphabet::getComplementCode( sequence[i] );	// and its reverse complement
	}
}
