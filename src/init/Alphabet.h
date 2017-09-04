#ifndef ALPHABET_H_
#define ALPHABET_H_

#include <cstring>	// e.g. std::strcmp
#include <iostream>	// e.g. std::cerr, std::cout, std::endl
#include <iomanip>  // e.g. std::setprecision

#include <ctype.h>	// e.g. tolower
#include <stdint.h>	// e.g. uint8_t
#include <stdlib.h>	// e.g. free

class Alphabet{

public:

	static void		init( char* alphabetType );
	static void		destruct();

	static size_t	getSize();
	static char* 	getAlphabet();
	static uint8_t	getCode( char base );				// get encoding for base
	static char		getBase( uint8_t code );			// get base from encoding
	static uint8_t	getComplementCode( uint8_t code );	// get complement encoding from encoding

private:

	static size_t	size_;					// alphabet size
	static char* 	alphabet_;				// alphabet bases ([A,C,G,T], [A,C,G,T,mC], ...)
	static char* 	complementAlphabet_;	// complementary alphabet bases ([T,G,C,A], [T,G,C,A,G], ...)
	static uint8_t* baseToCode_;			// convert base to encoding
	static char* 	codeToBase_;			// convert encoding to base
	static uint8_t* codeToComplementCode_;	// convert encoding to complement encoding
};

inline uint8_t Alphabet::getCode( char base ){
	return baseToCode_[( size_t )base];
}

inline char Alphabet::getBase( uint8_t code ){
	return codeToBase_[code];
}

inline uint8_t Alphabet::getComplementCode( uint8_t code ){
	return codeToComplementCode_[code];
}

#endif /* ALPHABET_H_ */
