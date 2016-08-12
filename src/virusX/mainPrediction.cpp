#include "GlobalPrediction.h"
#include "VirusHostInteractions.h"

int main( int nargs, char* args[] ){

	Global::init( nargs, args );

	VirusHostInteractions* interactions = new VirusHostInteractions( Global::inputDirectoryBaMMs, Global::extensionBaMMs );
	interactions->predict( Global::inputDirectorySeqs, Global::extensionSeqs );
	interactions->write( Global::outputDirectory );

	if( Global::verbose ){
		interactions->print();
	}

	Global::destruct();

	return 0;
}
