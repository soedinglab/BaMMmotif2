#ifndef MOTIFSET_H_
#define MOTIFSET_H_

#include "Motif.h"
#include "utils.h"

class MotifSet {

public:

	MotifSet();
	~MotifSet();

	std::vector<Motif*> getMotifs();		// get motifs
	size_t		 		getN();				// get number of motifs

	void 				print();			// print motifs to console
	void 				write();			// write motifs to files (basename.bmm)

private:

	std::vector<Motif*> motifs_;			// motifs
	size_t         		N_ = 0;				// number of motifs

};

#endif /* MOTIFSET_H_ */
