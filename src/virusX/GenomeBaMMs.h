#ifndef GENOMEBAMMS_H_
#define GENOMEBAMMS_H_

#include <list>

#include <dirent.h>	// e.g. opendir
#include <string.h>	// e.g. strcmp

#include "../shared/BackgroundModel.h"
#include "../shared/SequenceSet.h"
#include "../shared/utils.h"

#include "Global.h"

class GenomeBaMMs{

public:

	GenomeBaMMs();
	~GenomeBaMMs();

private:

	std::list<BackgroundModel*> bamms_;
};

#endif /* GENOMEBAMMS_H_ */
