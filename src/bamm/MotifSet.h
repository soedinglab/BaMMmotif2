#ifndef MOTIFSET_H_
#define MOTIFSET_H_

#include "Motif.h"
#include "../shared/utils.h"

class MotifSet {

public:

	MotifSet();
	~MotifSet();

	std::vector<Motif*> getMotifs();		// get motifs
	int			 		getN();				// get number of motifs

	void 				print();			// print motifs to console
	void 				write();			// write motifs to files (basename.bmm)

private:

	std::vector<Motif*> motifs_;			// motifs
	int          		N_ = 0;				// number of motifs
	std::vector<int>	Y_;					// contains 1 at position 0
											// and the number of oligomers y for increasing order k (from 0 to K_) at positions k+1
											// e.g. alphabet size_ = 4 and K_ = 2: Y_ = 1 4 16 64

};

#endif /* MOTIFSET_H_ */
