#ifndef SEQGENERATOR_H_
#define SEQGENERATOR_H_

#include "Alphabet.h"
#include "BackgroundModel.h"
#include "utils.h"
#include "Motif.h"
#include "Global.h"

class SeqGenerator {

	/*
	 * This class generates artificial sequences set for:
	 * 1. negative sequence set, by using k-mer frequencies;
	 * 2. simulated positive sequence set by inserting motif
	 *    into the negative sequence set
	 * 3. negative sequence set with motif patterns masked from the postive sequences
	 * Prerequisite:
	 * 1. sequences set for calculating k-mer frequencies
	 * 2. the order of k-mer for sampling
	 * 3. (optional) the motif model for generating pseudo-sequence set
	 */

public:
	SeqGenerator( std::vector<Sequence*> seqs,
					Motif* motif = NULL,
					size_t sOrder = 2 );
	~SeqGenerator();

	std::vector<std::unique_ptr<Sequence>> arti_negset( size_t fold );
	std::vector<std::unique_ptr<Sequence>> arti_posset_motif_embedded( size_t fold );
	std::vector<std::unique_ptr<Sequence>> arti_negset_motif_masked( float** r );


	void write( char* odir,
				std::string basename,
				size_t n,
				std::vector<std::unique_ptr<Sequence>> seqset );

private:

	void						calculate_kmer_frequency();
	std::unique_ptr<Sequence> 	negseq_dimer_freq( size_t L );
	std::unique_ptr<Sequence> 	posseq_motif_embedded(size_t L);
	std::unique_ptr<Sequence>	negseq_motif_masked( Sequence* posseq,
													size_t W,
													float* r );

	std::vector<Sequence*> 		seqs_;			// positive sequence set
	float**						freqs_;			// k-mer frequencies
	size_t** 					count_;			// k-mer counts
	Motif* 						motif_;			// the optimized motif
	size_t						sOrder_;		// the order of k-mers for
												// generating negative/pseudo
												// sequence set
	std::vector<size_t>			Y_;
};

#endif /* SEQGENERATOR_H_ */
