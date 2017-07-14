#ifndef MOTIFSET_H_
#define MOTIFSET_H_

#include "Motif.h"
#include "utils.h"

class MotifSet {

public:

	MotifSet( char* indir, size_t l_flank, size_t r_flank,
			size_t order, std::string tag );
	~MotifSet();

	std::vector<Motif*> getMotifs();		// get motifs
	size_t		 		getN();				// get number of motifs

	void 				print();			// print motifs to console
	void 				write( char* outdir ); // write motifs to file

private:

	std::vector<Motif*> motifs_;			// motifs
	size_t         		N_;					// number of motifs
	size_t         		l_flank_;			// size of the left flanking region
	size_t         		r_flank_;			// size of the right flanking region
	size_t         		K_;					// order of the motifs (for PWMs)
	char*				indir_;				// input directory for motif file
	std::string			tag_;				// indicates the motif format (PWMs,
											// binding sites or bamm patterns)
};

#endif /* MOTIFSET_H_ */
