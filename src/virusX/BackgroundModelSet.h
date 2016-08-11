#ifndef BACKGROUNDMODELSET_H_
#define BACKGROUNDMODELSET_H_

#include <fstream>	// e.g. std::ifstream
#include <list>

#include <dirent.h>		// e.g. opendir
#include <string.h>		// e.g. strcmp
#include <sys/stat.h>

#include "../shared/BackgroundModel.h"
#include "../shared/SequenceSet.h"
#include "../shared/utils.h"

#include "Global.h"

class BackgroundModelSet{

public:

	// learn background models from sequences
	BackgroundModelSet( char* inputDirectory, char* extension, int order, std::vector<float> alpha );
	// read in background models from files
	BackgroundModelSet( char* inputDirectory, char* extension );

	~BackgroundModelSet();

	std::list<BackgroundModel*>& getBackgroundModels();
	int getN();

	void print();
	void write( char* dir );

private:

	std::list<BackgroundModel*>	backgroundModels_;
	int          				N_ = 0;				// number of background models
};

#endif /* BACKGROUNDMODELSET_H_ */
