/*
 * Alphabet.cpp
 *
 *  Created on: Apr 12, 2016
 *      Author: wanwan
 */

#include <cstring>      // strlen

#include <stdlib.h>		// malloc, calloc, exit, free
#include <stdio.h>		// stderr
#include <ctype.h>		// tolower

#include "Alphabet.h"

void Alphabet::init( char const* alphabetType ){

	if( strcmp( alphabetType, "STANDARD" ) == 0 ){
		size_ = 4;
		alphabet_ = "ACGT";
		complementAlphabet_ = "TGCA";
	} else if( strcmp( alphabetType, "METHYLC" ) == 0 ){
		size_ = 5;
		alphabet_ = "ACGTM";
		complementAlphabet_ = "TGCAG";
	} else if( strcmp( alphabetType, "HYDROXYMETHYLC" ) == 0 ){
		size_ = 5;
		alphabet_ = "ACGTH";
		complementAlphabet_ = "TGCAG";
	} else if( strcmp( alphabetType, "EXTENDED" ) == 0 ){
		size_ = 6;
		alphabet_ = "ACGTMH";
		complementAlphabet_ = "TGCAGG";
	} else {
		fprintf( stderr, "Please provide one of the following alphabet types: "
		        "STANDARD, METHYLC, HYDROXYMETHYLC, EXTENDED.\n" );
		exit( -1 );
	}

	baseToCode_ = ( uint8_t* )calloc( 128, sizeof( uint8_t ) );
	for( unsigned int i = 0; i < size_; i++ ){
	  	baseToCode_[( int )alphabet_[i]] = i + 1;
	  	baseToCode_[( int )tolower( alphabet_[i] )] = i + 1;
	}
	codeToComplementCode_ = ( uint8_t* )calloc( size_+1, sizeof( uint8_t ) );
	for( unsigned int i = 0; i < size_; i++ ){
		codeToComplementCode_[i+1] = baseToCode_[( int )complementAlphabet_[i+1]];
	}
}

char const* Alphabet::getAlphabet(){
	return alphabet_;
}

char const* Alphabet::getComplementAlphabet(){
	return complementAlphabet_;
}

void Alphabet::destruct(){
	// free memory allocated to the pointers
	if( alphabet_ ) delete alphabet_;
	if( complementAlphabet_ ) delete complementAlphabet_;
	if( baseToCode_ ) free( baseToCode_ );
	if( codeToComplementCode_ ) free( codeToComplementCode_ );
}



