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
    enriched_kmer_expects_.resize(0);

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

        size_t seq_length = (Global::ss) ? seqSet[n]->getL() : ( seqSet[n]->getL()-1 ) / 2;

        for( size_t pos = 0; pos < seq_length; pos++ ) {

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
                encode = 0;
                step = 0;
                if( Global::verbose ){
                    std::cout << "\tNote: While counting k-mers, skipped undefined alphabet "
                              << Alphabet::getBase(base) << std::endl;
                }
            }

            if( step >= kmer_length_ ){
                size_t revcomp = encode2revcomp_[encode];
                size_t min_encode = ( encode <= revcomp ) ? encode : revcomp;
                //size_t min_encode = encode;
                // todo: there is a bug in double counting!!!
                assert( encode2index_[min_encode] < kmer_N_ );
                kmer_counts[encode2index_[min_encode]]++;
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
    enriched_kmer_expects_.resize( enriched_kmer_N_ );

    kmer_is_counted_ = true;
}

void KmerPredictor::scoreKmer(Motif* motif, BackgroundModel* bg) {

    assert( kmer_is_counted_ );

    size_t W = motif->getW();
    size_t K = motif->getK();

    /**
     * get log odds score between Vmotif and Vbg
     */
    motif->calculateLogS( bg->getV() );

    for( size_t idx = 0; idx < enriched_kmer_N_; idx++ ) {

        size_t encode = enriched_kmer_encodes_[idx];
        size_t revcomp = encode2revcomp_[encode];

        float best_score = 0.f;

        for (size_t i = 0; i < W + kmer_length_-1; i++ ){

            /**
             * calculate score for forward motif
             */
            // initialize y with a random (K+1)-mer
            size_t y_prev = 0;
            for (size_t j = 0; j < K+1; j++ ){
                y_prev += ( rand() % size2power_[1] ) * size2power_[K-j] ;
            }
            size_t y_prev2 = y_prev;

            float score = 0.f;
            //float score = 1.f;
            // calculate y for each position
            for (size_t pos = 0; pos < W; pos++) {

                if( i+pos >= W-1 and i+pos < kmer_length_+W-1 ){

                    size_t j = i+pos-W+1;
                    size_t base = ( encode % size2power_[j+1] ) / size2power_[j];
                    size_t y_new = y_prev % size2power_[K] * size2power_[1] + base;
                    y_prev = y_new;
                    score += motif->getS()[y_new][pos];
                    //score *= motif->getV()[K][y_new][pos];
                } else {
                    ;
                }
                if( i >= W-1 ){
                    size_t j = i-W+1;
                    size_t base = ( encode % size2power_[j+1] ) / size2power_[j];
                    y_prev = y_prev2 % size2power_[K] * size2power_[1] + base;
                    y_prev2 = y_prev;
                }
            }

            // take the best score within the motif region
            best_score = ( best_score > score) ? best_score : score;

            // sum up all the scores within the motif region
            //best_score += score;

            /**
             * calculate score for the reverse complement motif
             */
            y_prev = 0;
            for (size_t j = 0; j < K+1; j++ ){
                y_prev += ( rand() % size2power_[1] ) * size2power_[K-j] ;
            }
            y_prev2 = y_prev;

            score = 0.f;
            //score = 1.f;
            // calculate y for each position
            for (size_t pos = 0; pos < W; pos++) {
                if( i+pos >= W-1 and i+pos < kmer_length_+W-1 ){
                    size_t j = i+pos-W+1;
                    size_t base = ( revcomp % size2power_[j+1] ) / size2power_[j];
                    size_t y_new = y_prev % size2power_[K] * size2power_[1] + base;
                    y_prev = y_new;
                    score += motif->getS()[y_new][pos];
                    //score *= motif->getV()[K][y_new][pos];

                } else {
                    ;
                }
                if (i >= W-1){
                    size_t j = i-W+1;
                    size_t base = ( revcomp % size2power_[j+1] ) / size2power_[j];
                    y_prev = y_prev2 % size2power_[K] * size2power_[1] + base;
                    y_prev2 = y_prev;
                }
            }
            best_score = ( best_score > score) ? best_score : score;
            //best_score += score;
        }

        //float pred_score = expf( best_score ) / (1+expf( best_score ));
        //enriched_kmer_scores_[idx] = pred_score;
        enriched_kmer_scores_[idx] = best_score;
    }

}

void KmerPredictor::predictKmer( SequenceSet* seqSet,
                                 Motif *motif, BackgroundModel *bg) {

    assert( kmer_is_counted_ );

    size_t W = motif->getW();
    size_t K = motif->getK();

    // pre-calculate log odds scores given motif and bg model
    motif->calculateLogS( bg->getV() );
    float** s = motif->getS();
    std::vector<float> scores;
    size_t seqN = seqSet->getSequences().size();

    // convert all sequences to strings
    std::vector<std::string> seqs = seqSet->sequence2string();

    /**
     * store the log odds scores at all positions of each sequence
     */
    for( size_t n = 0; n < seqN; n++ ){

        size_t 	LW1 = seqSet->getSequences()[n]->getL() - W + 1;
        size_t* kmer = seqSet->getSequences()[n]->getKmer();
        float   score = 0.f;

        for( size_t i = 0; i < LW1; i++ ) {
            float logOdds = 0.0f;
            for (size_t j = 0; j < W; j++) {
                size_t y = kmer[i+j] % Global::A2powerK[K+1];
                logOdds += s[y][j];
            }
            // take all the log odds scores for MOPS model:
            score += expf(logOdds);
            //score += logOdds;
        }

        scores.push_back(score);
    }

    /**
     * Scan the whole sequences set and sum up scores
     * of which sequence contains motif
     */
    for( size_t idx = 0; idx < enriched_kmer_N_; idx++ ) {

        size_t      encode = enriched_kmer_encodes_[idx];
        std::string kmer_string = encode2string(encode);
        size_t      revcomp = encode2revcomp_[encode];
        std::string kmer_revcomp_string = encode2string(revcomp);

        /**
         * pre-calculate p_bg
         */
        float   p_bg = 1.f;
        size_t  k_bg = bg->getOrder();
        std::vector<size_t> kmers = string2kmer(kmer_string, k_bg);
        for (size_t j = 0; j < kmer_length_; j++) {
            size_t y = kmers[j];
            p_bg *= bg->getV()[k_bg][y];
        }

        float sum_score = 0.f;

        for (size_t n = 0; n < seqN; n++) {
            std::string sequence = seqs[n];

            if (sequence.find(kmer_string) != std::string::npos or
               (Global::ss and sequence.find(kmer_revcomp_string) != std::string::npos)) {
                sum_score += scores[n] * p_bg / (seqs[n].size() - W);
            }
        }
        enriched_kmer_expects_[idx] = sum_score;
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

std::vector<size_t> KmerPredictor::string2kmer(std::string kmer_string, size_t k_bg) {

    std::vector<size_t> kmer;
    kmer.resize(kmer_length_);

    for( size_t i = 0; i < kmer_length_; i++ ){
        kmer[i] = 0;
        //for( size_t k = i < Global::bgModelOrder ? i+1 : Global::bgModelOrder+1; k > 0; k-- ){
        for( size_t k = 0; k <= k_bg; k++ ){
            size_t base = ( i < k ) ?
                          ( size_t )rand() % Global::A2powerK[1] :
                          ( size_t )Alphabet::getCode(kmer_string[i-k])-1;
            kmer[i] += base * Global::A2powerK[k];
        }
    }

    return kmer;

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
    ofile << "kmer_IUPAC\tkmer_count\tkmer_score\texpected_kmer" << std::endl;
    for( size_t i = 0; i < enriched_kmer_N_; i++ ){
        ofile << encode2string(enriched_kmer_encodes_[i]) << '\t'
              << enriched_kmer_counts_[i] << '\t'
              << std::setprecision(4) << enriched_kmer_scores_[i] << '\t'
              << enriched_kmer_expects_[i] << std::endl;
    }
}