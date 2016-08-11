#include "BackgroundModelSet.h"
#include "Global.h"

#include "../shared/utils.h"

int main( int nargs, char* args[] ){

	Global::init( nargs, args );

//	BackgroundModelSet* bamms = new BackgroundModelSet( Global::inputDirectory, Global::extension, Global::modelOrder, Global::modelAlpha );
	BackgroundModelSet* bamms = new BackgroundModelSet( Global::inputDirectory, Global::extension );
	bamms->write( Global::outputDirectory );

	if( Global::verbose ){
		bamms->print();
	}

	Global::destruct();

	return 0;
}
