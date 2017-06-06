#ifndef SEQGENERATOR_H_
#define SEQGENERATOR_H_

#include "../shared/Alphabet.h"
#include "../shared/BackgroundModel.h"
#include "../shared/utils.h"
#include "Motif.h"
#include "ModelLearning.h"
#include "Global.h"

class SeqGenerator {

	/*
	 * This class generates artificial sequences set for:
	 * 1. negative sequence set, by using k-mer frequencies;
	 * 2. simulated positive sequence set by inserting motif
	 *    into the negative sequence set
	 */

public:
	SeqGenerator( std::vector<Sequence*> seqs, Motif* motif = NULL );
	~SeqGenerator();

	std::vector<std::unique_ptr<Sequence>> 	sample_negative_seqset( );
	std::vector<std::unique_ptr<Sequence>> 	sample_pseudo_seqset();
	std::vector<Sequence*>                  getSeqs();

	void									write( std::vector<std::unique_ptr<Sequence>> );

private:

	void									calculate_kmer_frequency();
	std::unique_ptr<Sequence> 				sample_negative_sequence( int L );
	std::unique_ptr<Sequence> 				sample_pseudo_sequence( int L );

	std::vector<Sequence*> 	seqs_;			// positive sequence set
	float**					freqs_;			// k-mer frequencies
	int** 					count_;			// k-mer counts
	Motif* 					motif_;			// the optimized motif

	std::vector<int>		Y_;				// contains 1 at position 0
											// and the number of oligomers y for increasing order k (from 0 to K_) at positions k+1
											// e.g. alphabet size_ = 4 and K_ = 2: Y_ = 1 4 16 64
};

#endif /* SEQGENERATOR_H_ */
