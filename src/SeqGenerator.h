#ifndef SEQGENERATOR_H_
#define SEQGENERATOR_H_

#include "Alphabet.h"
#include "BackgroundModel.h"
#include "utils.h"
#include "Motif.h"
#include "ModelLearning.h"
#include "Global.h"

class SeqGenerator {

	/*
	 * This class generates artificial sequences set for:
	 * 1. negative sequence set, by using k-mer frequencies;
	 * 2. simulated positive sequence set by inserting motif
	 *    into the negative sequence set
	 * Prerequisite:
	 * 1. sequences set for calculating k-mer frequencies
	 * 2. the order of k-mer for sampling
	 * 3. (optional) the motif model for generating pseudo-sequence set
	 */

public:
	SeqGenerator( std::vector<Sequence*> seqs, Motif* motif = NULL, size_t sOrder = Global::sOrder );
	~SeqGenerator();

	std::vector<std::unique_ptr<Sequence>> 	sample_negative_seqset( size_t fold );
	std::vector<std::unique_ptr<Sequence>> 	sample_pseudo_seqset( size_t fold );

	void				write_pseudoset( char* outdir, std::string basename );

private:

	void									calculate_kmer_frequency();
	std::unique_ptr<Sequence> 				sample_negative_sequence( size_t L );
	std::unique_ptr<Sequence> 				sample_pseudo_sequence( size_t L );

	std::vector<Sequence*> 	seqs_;			// positive sequence set
	float**					freqs_;			// k-mer frequencies
	size_t** 				count_;			// k-mer counts
	Motif* 					motif_;			// the optimized motif
	size_t					sOrder_;		// the order of k-mers for generating negative/pseudo sequence set
	std::vector<size_t>		Y_;				// contains 1 at position 0
											// and the number of oligomers y for increasing order k (from 0 to K_) at positions k+1
											// e.g. alphabet size_ = 4 and K_ = 2: Y_ = 1 4 16 64
};

#endif /* SEQGENERATOR_H_ */
