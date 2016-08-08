#include "GenomeBaMMs.h"
#include "Global.h"

int main( int nargs, char* args[] ){

	Global::init( nargs, args );

	GenomeBaMMs bamms;

	Global::destruct();

	return 0;
}
