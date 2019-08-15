//
// Created by wanwan on 13.08.19.
//

#include "KmerPredicter.h"

KmerPredicter::KmerPredicter( size_t kmer_length ){

    kmer_length_ = kmer_length;
    kmer_size_ = ipow( Alphabet::getSize(), kmer_length_ );
    enriched_kmer_scores_ = ( float* )calloc( kmer_size_, sizeof( float ) );
    enriched_kmer_ids_ = ( size_t* )calloc( kmer_size_, sizeof( size_t ) );
    enriched_kmer_counts_ = ( unsigned long* )calloc( kmer_size_, sizeof( unsigned long ) );
    size2power_ = ( size_t* )calloc( kmer_length_+1, sizeof( size_t ) );
    for (size_t i = 0; i < kmer_length_+1; i++){
        size2power_[i] = ipow( Alphabet::getSize(), i );
    }
    enriched_kmer_N_ = 0;
    kmer_is_counted_ = false;
}

KmerPredicter::~KmerPredicter() {
    free( enriched_kmer_scores_ );
    free( enriched_kmer_ids_ );
    free( enriched_kmer_counts_ );
    free( size2power_ );
}

void KmerPredicter::countKmer( std::vector<Sequence*> seqSet ) {

    std::vector<unsigned long> kmer_counts;
    kmer_counts.resize( kmer_size_ );
    for( size_t i = 0; i < kmer_size_; i++ ){
        kmer_counts[i] = 0;
    }

    for( size_t n = 0; n < seqSet.size(); n++ ){
        size_t id = 0;
        // size_t last_base = 0;
        // count the first k-mer from the starting base
        // for example: ATGCA <-> ID: 0*a^0 + 3*a^1 + 2*a^2 + 1*a^3 + 0*a^4
        for( size_t pos = 0; pos < kmer_length_; pos++ ){
            id += (seqSet[n]->getSequence()[pos]-1) * size2power_[pos];
        }
        kmer_counts[id]++;

        // count the k-mers till the end of the sequence
        for( size_t pos = kmer_length_; pos < seqSet[n]->getL(); pos++ ){
            id = id / size2power_[1] +
                    (seqSet[n]->getSequence()[pos]-1) * size2power_[kmer_length_-1];
            kmer_counts[id]++;
        }
    }

    /**
     * Only save k-mers with more than N occurrences
     */
    if( Global::kmerNCutoff ) {
        size_t n = 0;
        for (size_t idx = 0; idx < kmer_size_; idx++) {
            if (kmer_counts[idx] >= Global::kmerNCutoff ) {
                enriched_kmer_ids_[n] = idx;
                enriched_kmer_counts_[n] = kmer_counts[idx];
                n++;
            }
        }
        enriched_kmer_N_ = n;
    } else {
        enriched_kmer_N_ = kmer_size_;
    }
    kmer_is_counted_ = true;
}

void KmerPredicter::scoreKmer(Motif* motif, BackgroundModel* bg) {

    assert( kmer_is_counted_ );
    size_t W_ = motif->getW();
    size_t K_ = motif->getK();

    /**
     * get linear score between Vmotif and Vbg
     */
    motif->calculateLinearS( bg->getV() );

    if( kmer_length_ <= W_ ){
        /**
         * when k-mer is no longer than the motif
         */
        for( size_t idx = 0; idx < enriched_kmer_N_; idx++ ) {
            float best_score = 0.f;
            for (size_t i = 0; i <= W_ - kmer_length_; i++ ){
                float score = 1.f;
                // initialize y with a random (K+1)-mer
                size_t y_prev = 0;
                for (size_t j = 0; j < K_+1; j++ ){
                    y_prev += ( rand() % size2power_[1] ) * size2power_[K_-j] ;
                }

                // calculate y for each position
                for (size_t pos = 0; pos < kmer_length_; pos++) {
                    size_t base = ( enriched_kmer_ids_[idx] % size2power_[pos+1] ) / size2power_[pos];
                    size_t y_new = y_prev % size2power_[K_] * size2power_[1] + base;
                    score *= motif->getS()[y_new][pos+i];
                    y_prev = y_new;
                }
                best_score = ( best_score > score) ? best_score : score;
            }
            enriched_kmer_scores_[idx] = best_score;
        }
    } else {
        /**
         * when k-mer is longer than the motif
         */
        for( size_t idx = 0; idx < enriched_kmer_N_; idx++ ) {
            float best_score = 0.f;
            for (size_t i = 0; i < kmer_length_ - W_; i++ ){
                float score = 0.f;
                // initialize y with a random (K+1)-mer
                size_t y_prev = 0;
                for (size_t j = 0; j < K_+1; j++ ){
                    y_prev += ( rand() % size2power_[1] ) * size2power_[K_-j] ;
                }
                // calculate y for each position
                for (size_t pos = 0; pos < W_; pos++) {
                    size_t base = ( enriched_kmer_ids_[idx] % size2power_[pos+1] ) / size2power_[pos];
                    size_t y_new = y_prev % size2power_[K_] * size2power_[1] + base;
                    score *= motif->getS()[y_new][pos];
                    y_prev = y_new;
                }
                best_score = (best_score > score) ? best_score:score;
            }
            enriched_kmer_scores_[idx] = best_score;
        }
    }
}

std::string KmerPredicter::ID2String( size_t kmer_id ){

    std::string kmer_string;
    for( size_t pos = 0; pos < kmer_length_; pos++ ){
        size_t code = ( kmer_id % size2power_[pos+1] ) / size2power_[pos];
        kmer_string += Alphabet::getBase( code + 1 );
    }
    return kmer_string;
}

void KmerPredicter::writeKmerCounts(char *odir, std::string basename) {
    /**
     * save kmer counts in .kmercounts
     */
    std::string opath = std::string( odir )  + '/' + basename + ".kmercounts";
    std::ofstream ofile( opath );

    // add a header to the results
    ofile << "kmer_IUPAC\tkmer_count\tkmer_score" << std::endl;
    for( size_t i = 0; i < enriched_kmer_N_; i++ ){
        ofile << ID2String(enriched_kmer_ids_[i]) << '\t'
              << enriched_kmer_counts_[i] << '\t'
              << enriched_kmer_scores_[i] << std::endl;
    }
}