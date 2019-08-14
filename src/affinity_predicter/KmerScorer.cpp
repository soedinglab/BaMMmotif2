//
// Created by wanwan on 14.08.19.
//

#include "KmerScorer.h"

KmerScorer::KmerScorer( Motif* motif,
                        BackgroundModel* bg){

    motif_ = motif;
    bg_ = bg;
    W_ = motif_->getW();
    K_ = motif_->getK();
    kmer_length_ = Global::kmerLength;
    kmer_size_ = ipow( Alphabet::getSize(), kmer_length_ );
    size2power_ = new size_t[kmer_length_+1];
    for (size_t i = 0; i < kmer_length_+1; i++){
        size2power_[i] = ipow( Alphabet::getSize(), i );
    }

}

KmerScorer::~KmerScorer() {

    delete[] size2power_;

}

void KmerScorer::scoreKmer() {

    /**
     * get linear score between Vmotif and Vbg
     */
    motif_->calculateLinearS( bg_->getV() );


    if( kmer_length_ <= W_ ){
        /**
         * when kmer is no longer than the motif
         */
        for( size_t kmer_id = 0; kmer_id < kmer_size_; kmer_id++ ) {
            float best_score = 0.f;
            for (size_t i = 0; i < W_ - kmer_length_; i++ ){
                float score = 0.f;
                // initialize y with a random (K+1)-mer
                size_t y_prev = 0;
                for (size_t j = 0; j < K_+1; j++ ){
                    y_prev += size2power_[K_-j] * ( rand() % size2power_[1] );
                }
                // calculate y for each position
                for (size_t pos = 0; pos < kmer_length_; pos++) {
                    size_t base = kmer_id / size2power_[pos+1] % size2power_[1];
                    size_t y_new = y_prev % size2power_[K_] + base;
                    score *= motif_->getS()[y_new][pos+i];
                    y_prev = y_new;
                }
                best_score = (best_score > score) ? best_score:score;
            }
            kmer_scores_[kmer_id] = best_score;
        }
    } else {
        /**
         * when kmer is longer than the motif
         */
        for( size_t kmer_id = 0; kmer_id < kmer_size_; kmer_id++ ) {
            float best_score = 0.f;
            for (size_t i = 0; i < kmer_length_ - W_; i++ ){
                float score = 0.f;
                // initialize y with a random (K+1)-mer
                size_t y_prev = 0;
                for (size_t j = 0; j < K_+1; j++ ){
                    y_prev += size2power_[K_-j] * ( rand() % size2power_[1] );
                }
                // calculate y for each position
                for (size_t pos = 0; pos < W_; pos++) {
                    size_t base = kmer_id / size2power_[pos+1] % size2power_[1];
                    size_t y_new = y_prev % size2power_[K_] + base;
                    score *= motif_->getS()[y_new][pos];
                    y_prev = y_new;
                }
                best_score = (best_score > score) ? best_score:score;
            }
            kmer_scores_[kmer_id] = best_score;
        }
    }
}

float* KmerScorer::getKmerScores() {
    return kmer_scores_;
}

void KmerScorer::writeKmerScores(char *odir, std::string basename) {
    /**
     * save kmer scores in .kmerscores
     */
    std::string opath = std::string( odir )  + '/' + basename + ".kmerscores";
    std::ofstream ofile( opath );

    // add a header to the results
    ofile << "kmer_id\tkmer_count" << std::endl;
    for( size_t i = 0; i < kmer_size_; i++ ){
        ofile << i << '\t' << std::setprecision( 3 )
              << kmer_scores_[i] << std::endl;
    }
}