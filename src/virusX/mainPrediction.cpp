#include "GPrediction.h"
#include "../shared/BackgroundModelSetScore.h"

int main( int nargs, char* args[] ){

	Global::init( nargs, args );

	BackgroundModelSetScore* interactions = new BackgroundModelSetScore( Global::inputDirectoryBaMMs, Global::extensionBaMMs );
	interactions->predict( Global::inputDirectorySeqs, Global::extensionSeqs );

	interactions->write( Global::outputDirectory );

	if( Global::verbose ){
		interactions->print();
	}

	delete interactions;
	Global::destruct();

	return 0;
}
