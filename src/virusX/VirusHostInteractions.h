#ifndef VIRUSHOSTINTERACTIONS_H_
#define VIRUSHOSTINTERACTIONS_H_

#include <dirent.h>	// e.g. opendir

#include "../shared/BackgroundModelSet.h"
#include "../shared/SequenceSet.h"
#include "../shared/utils.h"

class VirusHostInteractions{

public:

	VirusHostInteractions( char* inputDirectoryBaMMs, char* extensionBaMMs );
	~VirusHostInteractions();

	void predict( char* inputDirectorySeqs, char* extensionSeqs );

	void print();
	void write( char* outputDirectory );

private:

	float score( BackgroundModel& bamm, SequenceSet& sequenceSet );

	BackgroundModelSet* bamms_;
	std::vector<std::vector<float>> posteriors_;
};

#endif /* VIRUSHOSTINTERACTIONS_H_ */
