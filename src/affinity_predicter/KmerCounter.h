//
// Created by wanwan on 13.08.19.
//

#ifndef SCOREKMER_H
#define SCOREKMER_H

#include "../Global/Global.h"

class KmerCounter {
    /**
     * This class is aimed for:
     * count k-mers from the given sequence set
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
    KmerCounter( std::vector<Sequence*> seqSet, size_t kmer_length );
    ~KmerCounter();

    // decode kmer ID back to kmer string
    std::string ID2String(const size_t kmer_id);
    void        countKmer();
    size_t*     getKmerCounts();
    void        writeKmerCounts( char* odir, std::string basename );

private:

    std::vector<Sequence*> seqSet_;
    size_t  kmer_length_;
    size_t  kmer_size_;
    size_t* kmer_counts_;
    size_t* size2power_;

};


#endif //SCOREKMER_H
