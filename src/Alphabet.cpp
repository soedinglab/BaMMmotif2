/*
 * Alphabet.cpp
 *
 *  Created on: Apr 12, 2016
 *      Author: wanwan
 */

#include <cstring>      /* strlen */

#include "Alphabet.h"

void Alphabet::init( char const* alphabetString ){

    size_ = strlen( alphabetString );
    alphabet_ = ( char* )malloc( sizeof( char ) );
}



