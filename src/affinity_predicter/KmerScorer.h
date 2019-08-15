//
// Created by wanwan on 14.08.19.
//

#ifndef KMERSCORER_H
#define KMERSCORER_H

#include "../Global/Global.h"
#include "../init/Motif.h"
#include "../init/BackgroundModel.h"

class KmerScorer {

    /**
     * This class is aimed for:
     * score k-mers given model and sequence set
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
    KmerScorer( Motif* motif,
                BackgroundModel* bg);
    ~KmerScorer();

    void        scoreKmer();
    float*      getKmerScores();
    void        writeKmerScores( char* odir, std::string basename );

private:

    Motif*  motif_;
    BackgroundModel* bg_;
    size_t  W_;
    size_t  K_;
    size_t  kmer_length_;
    size_t  kmer_size_;
    float*  kmer_scores_;
    size_t* size2power_;

};


#endif //KMERSCORER_H
