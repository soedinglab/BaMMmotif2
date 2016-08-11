#include "Alphabet.h"

int 		Alphabet::size_ = 0; 				// alphabet size default that can be checked for
char*		Alphabet::alphabet_;				// alphabet bases ([A,C,G,T], [A,C,G,T,mC], ...)
char*		Alphabet::complementAlphabet_;		// complementary alphabet bases ([T,G,C,A], [T,G,C,A,G], ...)
uint8_t*	Alphabet::baseToCode_;				// convert base to encoding
char*		Alphabet::codeToBase_;				// convert encoding to base
uint8_t*	Alphabet::codeToComplementCode_;	// convert encoding to complement encoding

void Alphabet::init( char* alphabetType ){

	if( strcmp( alphabetType, "STANDARD" ) == 0 ){
		size_ = 4;
		alphabet_ = strdup( "ACGT" );
		complementAlphabet_ = strdup( "TGCA" );
	} else if( strcmp( alphabetType, "METHYLC" ) == 0 ){
		size_ = 5;
		alphabet_ = strdup( "ACGTM" );
		complementAlphabet_ = strdup( "TGCAG" );
	} else if( strcmp( alphabetType, "HYDROXYMETHYLC" ) == 0 ){
		size_ = 5;
		alphabet_ = strdup( "ACGTH" );
		complementAlphabet_ = strdup( "TGCAG" );
	} else if( strcmp( alphabetType, "EXTENDED" ) == 0 ){
		size_ = 6;
		alphabet_ = strdup( "ACGTMH" );
		complementAlphabet_ = strdup( "TGCAGG" );
	} else {
		std::cerr << "Error: Correct alphabet type to STANDARD, METHYLC, HYDROXYMETHYLC, or EXTENDED" << std::endl;
		exit( -1 );
	}

	baseToCode_ = ( uint8_t* )calloc( 128, sizeof( uint8_t ) );
	codeToBase_ = ( char* )calloc( 128, sizeof( char ) );
	for( int i = 0; i < size_; i++ ){
	  	// base:		N	A	C	G	T	...
	  	// encoding:	0	1	2	3	4	...
	  	baseToCode_[( int )alphabet_[i]] = ( uint8_t )( i + 1 );
	  	baseToCode_[( int )tolower( alphabet_[i] )] = ( uint8_t )( i + 1 );		// for lower case
	  	codeToBase_[i+1] = alphabet_[i];
	}

	codeToComplementCode_ = ( uint8_t* )calloc( size_+1, sizeof( uint8_t ) );
	for( int i = 0; i < size_; i++ )
		codeToComplementCode_[i+1] = baseToCode_[( int )complementAlphabet_[i]];
}

void Alphabet::destruct(){

	if( alphabet_ ){
		free( alphabet_ );
	}
	if( complementAlphabet_ ){
		free( complementAlphabet_ );
	}
	if( baseToCode_ ){
		free( baseToCode_ );
	}
	if( codeToBase_ ){
		free( codeToBase_ );
	}
	if( codeToComplementCode_ ){
		free( codeToComplementCode_ );
	}
}

char* Alphabet::getAlphabet(){
	return alphabet_;
}

char* Alphabet::getComplementAlphabet(){
	return complementAlphabet_;
}

int Alphabet::getSize(){
	return size_;
}

void Alphabet::setSize( int size ){
	size_ = size;
}
