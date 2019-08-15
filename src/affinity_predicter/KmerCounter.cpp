//
// Created by wanwan on 13.08.19.
//

#include "KmerCounter.h"

KmerCounter::KmerCounter( std::vector<Sequence*> seqSet ){

    seqSet_ = seqSet;
    kmer_length_ = Global::kmerLength;
    kmer_size_ = ipow( Alphabet::getSize(), kmer_length_ );
    kmer_counts_ = ( size_t* )calloc( kmer_size_, sizeof( size_t ) );
    for( size_t i = 0; i < kmer_size_; i++ ){
        kmer_counts_[i] = 0;
    }
    size2power_ = ( size_t* )calloc( kmer_length_+1, sizeof( size_t ) );
    for (size_t i = 0; i < kmer_length_+1; i++){
        size2power_[i] = ipow( Alphabet::getSize(), i );
    }

}

KmerCounter::~KmerCounter() {
    free( kmer_counts_ );
    free( size2power_ );
}


void KmerCounter::countKmer() {

    for( size_t n = 0; n < seqSet_.size(); n++ ){
        size_t id = 0;
        // size_t last_base = 0;
        // count the first k-mer from the starting base
        // for example: ATGCA <-> ID: 0*a^0 + 3*a^1 + 2*a^2 + 1*a^3 + 0*a^4
        for( size_t pos = 0; pos < kmer_length_; pos++ ){
            id += (seqSet_[n]->getSequence()[pos]-1) * size2power_[pos];
        }
        kmer_counts_[id]++;

        // count the k-mers till the end of the sequence
        for( size_t pos = kmer_length_; pos < seqSet_[n]->getL(); pos++ ){
            id = id / size2power_[1] +
                    (seqSet_[n]->getSequence()[pos]-1) * size2power_[kmer_length_-1];
            kmer_counts_[id]++;
        }
    }

}

std::string KmerCounter::ID2String( size_t kmer_id ){

    std::string kmer_string;
    for( size_t pos = 0; pos < kmer_length_; pos++ ){
        size_t code = ( kmer_id % size2power_[pos+1] ) / size2power_[pos];
        kmer_string += Alphabet::getBase( code + 1 );
    }
    return kmer_string;
}

size_t* KmerCounter::getKmerCounts() {
    return kmer_counts_;
}

void KmerCounter::writeKmerCounts(char *odir, std::string basename) {

    /**
     * save kmer counts in .kmercounts
     */
    std::string opath = std::string( odir )  + '/' + basename + ".kmercounts";
    std::ofstream ofile( opath );

    // add a header to the results
    ofile << "kmer_id\tkmer_IUPAC\tkmer_count" << std::endl;
    for( size_t i = 0; i < kmer_size_; i++ ){
        ofile << i << '\t'
              << ID2String(i) << '\t'
              << kmer_counts_[i] << std::endl;
    }
}