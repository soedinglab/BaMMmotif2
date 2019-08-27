//
// Created by wanwan on 13.08.19.
//

#ifndef SCOREKMER_H
#define SCOREKMER_H

#include "../Global/Global.h"
#include "../init/Motif.h"
#include "../init/BackgroundModel.h"

class KmerPredictor {

    /**
     * This class is aimed for:
     * count k-mers from the given sequence set
     *
     * k-mers are encoded with IDs
     * Alphabet: ACGT...
     * base     A   C   G   T
     * code     0   1   2   3
     * Alphabet size a (here a=4)
     * for example: ATGCA <-> ID: 0*a^0 + 3*a^1 + 2*a^2 + 1*a^3 + 0*a^4
     */

public:

    // define de-/constructor
    KmerPredictor( size_t kmer_length );
    ~KmerPredictor();

    // decode k-mer ID back to k-mer string
    std::string         encode2string(unsigned long kmer_encode);
    void                countKmer( std::vector<Sequence*> seqSet );
    void                scoreKmer(Motif* motif, BackgroundModel* bg);
    void                writeKmerStats(char *odir, std::string basename);
    void                calcRevComp();

private:

    size_t  kmer_length_;
    size_t  kmer_size_;
    size_t* size2power_;
    unsigned long* kmer2encode_;
    unsigned long* encode2revcomp_;
    unsigned long* encode2index_;
    size_t  kmer_N_;
    bool    kmer_is_counted_;

    std::vector<unsigned long> index2encode_;

    std::vector<float>  enriched_kmer_scores_;
    std::vector<size_t> enriched_kmer_encodes_;
    std::vector<size_t> enriched_kmer_counts_;

    size_t  enriched_kmer_N_;




};


#endif //SCOREKMER_H
