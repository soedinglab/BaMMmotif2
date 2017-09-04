#ifndef MOTIFSET_H_
#define MOTIFSET_H_

#include "Motif.h"
#include "../refinement/utils.h"

class MotifSet {

public:

	MotifSet( char* indir, size_t l_flank, size_t r_flank, std::string tag,
			  SequenceSet* posSet = NULL, float* f_bg = NULL, size_t order = 2, std::vector<float>  = {} );
	~MotifSet();

	std::vector<Motif*> getMotifs();		// get motifs
	size_t		 		getN();				// get number of motifs

	void 				print();			// print motifs to console
	void 				write( char* outdir ); // write motifs to file

private:

	std::vector<Motif*> motifs_;			// motifs
	size_t         		N_;					// number of motifs

};

#endif /* MOTIFSET_H_ */
