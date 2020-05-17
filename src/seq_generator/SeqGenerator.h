#ifndef SEQGENERATOR_H_
#define SEQGENERATOR_H_

#include <random>
#include "../init/Alphabet.h"
#include "../init/BackgroundModel.h"
#include "../init/Motif.h"
#include "../Global/Global.h"

class SeqGenerator {

	/**
	 * This class generates artificial sequences sets as listed:
	 * 1. negative sequence set, by using s-mer frequencies;
	 * 2. simulated positive sequence set, by inserting motif into the negative sequence set;
	 * 3. negative sequence set with motif patterns masked from the positive sequences.
	 * Prerequisites:
	 * 1. sequences set for calculating s-mer frequencies;
	 * 2. the order of s-mer for sampling, s = 2 by default;
	 * 3. (optional) the motif model for embedding or masking.
	 */

public:
	SeqGenerator( std::vector<Sequence*> seqs, Motif* motif = NULL );

	~SeqGenerator();

	std::vector<std::unique_ptr<Sequence>> sample_bgset_by_fold(size_t fold);
    std::vector<std::unique_ptr<Sequence>> sample_bgset_by_num(size_t negN, size_t maxL);
	std::vector<std::unique_ptr<Sequence>> arti_posset_motif_embedded(size_t at);
	std::vector<std::unique_ptr<Sequence>> seqset_with_motif_masked(float **r);
    std::vector<std::unique_ptr<Sequence>> seqset_with_best_motif_masked(float **r);

	void write( char* odir,
				std::string basename,
				std::vector<std::unique_ptr<Sequence>> seqset );

    void write( char* odir,
                std::string basename,
                std::vector<Sequence*> seqset );
private:

	void						calculate_kmer_frequency();
    void                        rescale_kmer_frequency( Sequence* refSeq );

	std::unique_ptr<Sequence> 	bg_sequence( size_t L, size_t count=0);
    std::unique_ptr<Sequence>   bgseq_on_rescaled_v( Sequence* refSeq, size_t count=0 );
    std::unique_ptr<Sequence> 	raw_seq(Sequence *refSeq);
	std::unique_ptr<Sequence> 	posseq_motif_embedded( Sequence* seq, size_t at );
	std::unique_ptr<Sequence>	posseq_motif_masked(Sequence* seq, size_t W, float *r);
    std::unique_ptr<Sequence>	posseq_best_motif_masked(Sequence* seq, size_t n, size_t W, float *r);

	std::vector<Sequence*> 		seqs_;			// positive sequence set

    float**                     v_;             // k-mer conditional probabilities
	size_t** 					n_;			    // k-mer counts
    float**                     range_bar_;     // store cumulated sum of k-mers

    // for re-scaling the background model
    float**                     v_seq_;         // sequence-specific k-mer conditional probabilities
    size_t** 					n_seq_;			// sequence-specific k-mer counts

    float*                      A_;             // pseudo-parameter for k-mer counting
	Motif* 						motif_;			// the optimized motif
	size_t						sOrder_;	    // the order of k-mers for generating negative/pseudo sequence set
    float                       q_;             // portion of sequences in the set that are masked/embedded with the motif
    bool                        genericNeg_;    // flag for generating sequence specific negative sequences

    size_t                      N_;             // input sequence number
    bool                        kmer_freq_is_calculated_;
    bool                        kmer_freq_is_rescaled_;

};


#endif /* SEQGENERATOR_H_ */
