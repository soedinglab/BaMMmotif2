//
// Created by wanwan on 13.08.19.
//

#include "ScoreKmer.h"

ScoreKmer::ScoreKmer( std::vector<Sequence*> seqSet, size_t kmer_length ){

    seqSet_ = seqSet;
    kmer_length_ = kmer_length;
    kmer_size_ = ipow( Alphabet::getSize(), kmer_length_ );
    kmer_counts_ = new size_t[kmer_size_];
    for( size_t i = 0; i < kmer_size_; i++ ){
        kmer_counts_[i] = 0;
    }
    size2power_ = new size_t[kmer_length_+1];
    for (size_t i = 0; i < kmer_length_+1; i++){
        size2power_[i] = ipow( Alphabet::getSize(), i );
    }

}

ScoreKmer::~ScoreKmer() {
    delete[] kmer_counts_;
    delete[] size2power_;
    delete[] kmer_ids_;
}


void  ScoreKmer::CountKmer() {
    for( size_t n = 0; n < seqSet_.size(); n++ ){

        size_t id = 0;
        size_t last_base = 0;

        // count the first k-mer from the starting base
        for( size_t pos = 0; pos < kmer_length_; pos++ ){
            id += seqSet_[n]->getSequence()[pos] * size2power_[pos+1];
        }
        kmer_counts_[id]++;
        last_base = seqSet_[n]->getSequence()[kmer_length_-1];

        // count the k-mers till the end of the sequence
        for( size_t pos = kmer_length_; pos < seqSet_[n]->getL(); pos++ ){
            id = id / 4 - last_base +  seqSet_[n]->getSequence()[pos] * size2power_[kmer_length_];
            kmer_counts_[id]++;
        }

    }
}

std::string ScoreKmer::ID2String(const size_t kmer_id){

    std::string kmer_string;
    for( size_t pos = 0; pos < kmer_length_; pos++ ){
        size_t code = kmer_id / size2power_[pos+1] % Alphabet::getSize();
        kmer_string += Alphabet::getBase( code+1 );
    }

    return kmer_string;
}
