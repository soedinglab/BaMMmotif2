#include "Global.h"

int main( int nargs, char* args[] ){

	Global::init( nargs, args );

//	BackgroundModel* bgModel = new BackgroundModel();

	Global::destruct();

	return 0;
}
