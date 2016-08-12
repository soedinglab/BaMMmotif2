#include "GlobalBuilding.h"
#include "../shared/BackgroundModelSet.h"

int main( int nargs, char* args[] ){

	Global::init( nargs, args );

	BackgroundModelSet* bamms = new BackgroundModelSet( Global::inputDirectory, Global::extension, Global::modelOrder, Global::modelAlpha );
	bamms->write( Global::outputDirectory );

	if( Global::verbose ){
		bamms->print();
	}

	Global::destruct();

	return 0;
}
