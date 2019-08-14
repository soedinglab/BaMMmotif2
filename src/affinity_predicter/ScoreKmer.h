//
// Created by wanwan on 13.08.19.
//

#ifndef SCOREKMER_H
#define SCOREKMER_H


#include "../init/Motif.h"
#include "../init/BackgroundModel.h"
#include "../Global/Global.h"

class ScoreKmer {
    /**
     * This class is aimed for:
     * 1) count k-mers from the given sequence set
     * 2) score k-mers given model and sequence set
     *
     * k-mers are encoded with IDs
     * Alphabet: ACGT...
     * base     A   C   G   T
     * code     0   1   2   3
     * Alphabet size a (here a=4)
     * for example: ATGCA <-> ID: 0*a^1 + 3*a^2 + 2*a^3 + 1*a^4 + 0*a^5
     */

public:

    // define de-/constructor
    ScoreKmer( std::vector<Sequence*> seqSet, size_t kmer_length );
    ~ScoreKmer();

    // decode kmer ID back to kmer string
    std::string ID2String(const size_t kmer_id);
    void CountKmer();


private:

    std::vector<Sequence*> seqSet_;
    size_t kmer_length_;
    size_t kmer_size_;
    size_t* kmer_ids_;
    size_t* kmer_counts_;
    size_t* size2power_;

};


#endif //SCOREKMER_H
