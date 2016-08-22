#include "GPrediction.h"
#include "../shared/BackgroundModelSetScore.h"

int main( int nargs, char* args[] ){

	Global::init( nargs, args );

	BackgroundModelSetScore* interactions = new BackgroundModelSetScore( Global::inputDirectoryBaMMs, Global::extensionBaMMs );
	interactions->score( Global::inputDirectorySeqs, Global::extensionSeqs, Global::outputDirectory );

	delete interactions;
	Global::destruct();

	return 0;
}
