//
// Created by wanwan on 13.08.19.
//

#ifndef SCOREKMER_H
#define SCOREKMER_H

#include "../Global/Global.h"
#include "../init/Motif.h"
#include "../init/BackgroundModel.h"

class KmerPredicter {
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
    KmerPredicter( size_t kmer_length );
    ~KmerPredicter();

    // decode kmer ID back to kmer string
    std::string         ID2String( size_t kmer_id );
    void                countKmer( std::vector<Sequence*> seqSet );
    void                scoreKmer(Motif* motif, BackgroundModel* bg);
    void                writeKmerCounts( char* odir, std::string basename );

private:

    size_t  kmer_length_;
    size_t  kmer_size_;
    float*  enriched_kmer_scores_;
    size_t* enriched_kmer_ids_;
    unsigned long* enriched_kmer_counts_;
    size_t  enriched_kmer_N_;
    size_t* size2power_;
    bool    kmer_is_counted_;

};


#endif //SCOREKMER_H
