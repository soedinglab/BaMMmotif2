#ifndef GENOMEBAMMS_H_
#define GENOMEBAMMS_H_

#include <list>

#include <dirent.h>		// e.g. opendir
#include <string.h>		// e.g. strcmp
#include <sys/stat.h>

#include "../shared/BackgroundModel.h"
#include "../shared/SequenceSet.h"
#include "../shared/utils.h"

#include "Global.h"

class GenomeBaMMs{

public:

	GenomeBaMMs();
	~GenomeBaMMs();

	void print();

private:

	std::list<BackgroundModel*> bamms_;
};

#endif /* GENOMEBAMMS_H_ */
