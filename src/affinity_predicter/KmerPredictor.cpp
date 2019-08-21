//
// Created by wanwan on 13.08.19.
//

#include "KmerPredictor.h"

KmerPredictor::KmerPredictor( size_t kmer_length ){

    kmer_length_ = kmer_length;
    kmer_size_ = 1;
    for (size_t i = 0; i < kmer_length_; i++){
        kmer_size_ += ( Alphabet::getSize() - 1 ) * ipow( Alphabet::getSize(), i );
    }

    kmer2encode_ = ( unsigned long* )calloc( kmer_size_, sizeof( unsigned long ) );
    encode2revcomp_ = ( unsigned long* )calloc( kmer_size_, sizeof( unsigned long ) );
    encode2index_ = ( unsigned long* )calloc( kmer_size_, sizeof( unsigned long ) );

    size2power_ = ( size_t* )calloc( kmer_length_+1, sizeof( size_t ) );
    for (size_t i = 0; i < kmer_length_+1; i++){
        size2power_[i] = ipow( Alphabet::getSize(), i );
    }
    kmer_N_ = 0;
    enriched_kmer_N_ = 0;
    index2encode_.resize(0);

    enriched_kmer_encodes_.resize(0);
    enriched_kmer_counts_.resize(0);
    enriched_kmer_scores_.resize(0);
    kmer_is_counted_ = false;
}

KmerPredictor::~KmerPredictor() {

    free( kmer2encode_ );
    free( encode2revcomp_ );
    free( encode2index_ );
    free( size2power_ );
}

void KmerPredictor::countKmer( std::vector<Sequence*> seqSet ) {

    /**
     * encode reverse complement of each k-mer
     */
    calcRevComp();

    std::vector<size_t> kmer_counts;
    kmer_counts.resize( kmer_N_ );
    for( size_t i = 0; i < kmer_N_; i++ ){
        kmer_counts[i] = 0;
    }

    for( size_t n = 0; n < seqSet.size(); n++ ){

        size_t encode = 0;
        size_t step = 0;

        for( size_t pos = 0; pos < seqSet[n]->getL(); pos++ ) {

            uint8_t base = seqSet[n]->getSequence()[pos];
            // to avoid exceeding
            base = Alphabet::getCode(Alphabet::getBase( base ));

            if ( base > 0 ) {
                if( step < kmer_length_ ){
                    encode += (base-1) * size2power_[step];
                } else {
                    encode = encode / size2power_[1] +
                            (base-1) * size2power_[kmer_length_-1];
                }
                step++;
            } else {
                //std::cout << "\tNote: Skipped undefined alphabet "
                //          << Alphabet::getBase(base) << std::endl;
                encode = 0;
                step = 0;
            }

            if( step >= kmer_length_ ){
                size_t new_encode;
                size_t revcomp = encode2revcomp_[encode];
                new_encode = ( encode <= revcomp ) ? encode : revcomp;
                assert( encode2index_[new_encode] < kmer_N_ );
                kmer_counts[encode2index_[new_encode]]++;
            }
        }
    }

    /**
     * Only save k-mers with more than N_cutoff occurrences
     */
    size_t k = 0;
    for (size_t index = 0; index < kmer_N_; index++) {
        if (kmer_counts[index] >= Global::kmerNCutoff ) {
            enriched_kmer_encodes_.push_back( index2encode_[index] );
            enriched_kmer_counts_.push_back( kmer_counts[index] );
            k++;
        }
    }
    enriched_kmer_N_ = k;
    enriched_kmer_scores_.resize( enriched_kmer_N_ );

    kmer_is_counted_ = true;
}

void KmerPredictor::scoreKmer(Motif* motif, BackgroundModel* bg) {

    assert( kmer_is_counted_ );

    size_t W_ = motif->getW();
    size_t K_ = motif->getK();

    /**
     * get linear score between Vmotif and Vbg
     */
    motif->calculateLogS( bg->getV() );

    /**
     * when k-mer is no longer than the motif
     */
    for( size_t idx = 0; idx < enriched_kmer_N_; idx++ ) {

        //std::cout << encode2string(enriched_kmer_encodes_[idx]) << ": ";
        size_t encode = enriched_kmer_encodes_[idx];
        size_t revcomp = encode2revcomp_[encode];

        float best_score = 0.f;

        for (size_t i = 0; i < W_ + kmer_length_-1; i++ ){

            float score = 0.f;
            // initialize y with a random (K+1)-mer
            size_t y_prev = 0;
            for (size_t j = 0; j < K_+1; j++ ){
                y_prev += ( rand() % size2power_[1] ) * size2power_[K_-j] ;
            }
            size_t y_prev2 = y_prev;

            // calculate y for each position
            for (size_t pos = 0; pos < W_; pos++) {
                if( i+pos >= W_-1 and i+pos < kmer_length_+W_-1 ){
                    size_t j = i+pos-W_+1;
                    size_t base = ( encode % size2power_[j+1] ) / size2power_[j];
                    size_t y_new = y_prev % size2power_[K_] * size2power_[1] + base;
                    y_prev = y_new;
                    score += motif->getS()[y_new][pos];
                } else {

                }
                if (i >= W_-1){
                    size_t j = i-W_+1;
                    size_t base = ( encode % size2power_[j+1] ) / size2power_[j];
                    y_prev = y_prev2 % size2power_[K_] * size2power_[1] + base;
                    y_prev2 = y_prev;
                }
            }
            best_score = ( best_score > score) ? best_score : score;

            // calculate score for the reverse complement
            score = 0.f;
            y_prev = 0;
            for (size_t j = 0; j < K_+1; j++ ){
                y_prev += ( rand() % size2power_[1] ) * size2power_[K_-j] ;
            }
            y_prev2 = y_prev;
            // calculate y for each position
            for (size_t pos = 0; pos < W_; pos++) {
                if( i+pos >= W_-1 and i+pos < kmer_length_+W_-1 ){
                    size_t j = i+pos-W_+1;
                    size_t base = ( revcomp % size2power_[j+1] ) / size2power_[j];
                    size_t y_new = y_prev % size2power_[K_] * size2power_[1] + base;
                    y_prev = y_new;
                    score += motif->getS()[y_new][pos];
                } else {

                }
                if (i >= W_-1){
                    size_t j = i-W_+1;
                    size_t base = ( revcomp % size2power_[j+1] ) / size2power_[j];
                    y_prev = y_prev2 % size2power_[K_] * size2power_[1] + base;
                    y_prev2 = y_prev;
                }
            }
            best_score = ( best_score > score) ? best_score : score;
        }

        enriched_kmer_scores_[idx] = best_score;
    }

}

std::string KmerPredictor::encode2string(unsigned long kmer_encode){

    std::string kmer_string;
    for( size_t pos = 0; pos < kmer_length_; pos++ ){
        size_t code = ( kmer_encode % size2power_[pos+1] ) / size2power_[pos];
        kmer_string += Alphabet::getBase( code + 1 );
    }
    return kmer_string;
}

void KmerPredictor::calcRevComp() {

    size_t index = 0;
    index2encode_.resize(0);

    for (size_t i = 0; i < kmer_size_; i++ ){
        size_t revcomp = 0;
        for (size_t j = 0; j < kmer_length_; j++ ){
            revcomp += ( size2power_[1] - 1 - (i % size2power_[j+1]) / size2power_[j] )
                       * size2power_[kmer_length_-j-1];
        }

        encode2revcomp_[i] = revcomp;

        if( i <= revcomp ){
            index2encode_.push_back(i);
            encode2index_[i] = index;
            kmer2encode_[index] = i;
            index ++;
        }
    }

    kmer_N_ = index;

    if( Global::verbose ) {
        /**
         * Printout for checking
         */
        std::cout << "There are " << kmer_N_ << " kmer indices.\n";
        for (size_t index = 0; index < kmer_N_; index++) {
            size_t encode = index2encode_[index];
            size_t revcomp = encode2revcomp_[encode];
            std::cout << encode2string(encode) << '\t'
                      << encode2string(revcomp) << '\t'
                      << index << '\t' << encode << '\t' << revcomp << '\t'
                      << encode2index_[encode] << '\t'
                      << encode2string(kmer2encode_[index]) << '\t'
                      << std::endl;
        }
    }


}

void KmerPredictor::writeKmerStats(char *odir, std::string basename) {

    /**
     * save kmer counts in .kmerstats
     */
    std::string opath = std::string( odir )  + '/' + basename + ".kmerstats";
    std::ofstream ofile( opath );

    // add a header to the results
    ofile << "kmer_IUPAC\tkmer_count\tkmer_score" << std::endl;
    for( size_t i = 0; i < enriched_kmer_N_; i++ ){
        ofile << encode2string(enriched_kmer_encodes_[i]) << '\t'
              << enriched_kmer_counts_[i] << '\t'
              << enriched_kmer_scores_[i] << std::endl;
    }
}