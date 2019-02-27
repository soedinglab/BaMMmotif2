#include "SeqGenerator.h"

SeqGenerator::SeqGenerator( std::vector<Sequence*> seqs, Motif* motif ){

	seqs_       = seqs;
	sOrder_     = Global::sOrder;
    motif_      = motif;
    q_          = Global::q;
    genericNeg_ = Global::genericNeg;
    
    v_ = ( float** )calloc( sOrder_+1, sizeof( float* ) );
	n_ = ( size_t** )calloc( sOrder_+1, sizeof( size_t* ) );
    range_bar_ = ( float** )calloc( sOrder_+1, sizeof( float* ) );
    v_seq_ = ( float** )calloc( sOrder_+1, sizeof( float* ) );
    n_seq_ = ( size_t** )calloc( sOrder_+1, sizeof( size_t* ) );
	for( size_t k = 0; k < sOrder_+1; k++ ) {
        v_[k] = (float *) calloc(Global::A2powerK[k + 1], sizeof(float));
        n_[k] = (size_t *) calloc(Global::A2powerK[k + 1], sizeof(size_t));
        range_bar_[k] = (float *) calloc(Global::A2powerK[k + 1], sizeof(float));
        v_seq_[k] = (float *) calloc(Global::A2powerK[k + 1], sizeof(float));
        n_seq_[k] = (size_t *) calloc(Global::A2powerK[k + 1], sizeof(size_t));
    }

    A_ = ( float* )calloc( sOrder_+1, sizeof( float ) );
    for( size_t k = 0; k < sOrder_+1; k++ ){
        A_[k] = 20.f;
    }

    rngx_.seed( 42 );
    srand( 42 );

    kmer_freq_is_calculated_    = false;
    kmer_freq_is_rescaled_      = false;

    N_ = seqs_.size();

}

SeqGenerator::~SeqGenerator(){

	for( size_t k = 0; k < sOrder_+1; k++ ){
        free( v_[k] );
		free( n_[k] );
        free( range_bar_[k] );
        free( v_seq_[k] );
        free( n_seq_[k] );
	}

    free( v_ );
	free( n_ );
    free( range_bar_ );
    free( v_seq_ );
	free( n_seq_ );

    free( A_ );

}

void SeqGenerator::calculate_kmer_frequency(){

	// reset counts for k-mers
	for( size_t k = 0; k < sOrder_+1; k++ ){
		for( size_t y = 0; y < Global::A2powerK[k+1]; y++ ){
			n_[k][y] = 0;
		}
	}

	// count k-mers
	for( size_t i = 0; i < seqs_.size(); i++ ){
		size_t L = seqs_[i]->getL();
		size_t* kmer = seqs_[i]->getKmer();
		for( size_t k = 0; k < sOrder_+1; k++ ){
			for( size_t j = k; j < L; j++ ){
				// extract (k+1)-mer
				size_t y = kmer[j] % Global::A2powerK[k+1];
				// count (k+1)mer
				n_[k][y]++;
			}
		}
	}

	// calculate frequencies
	size_t normFactor = 0;
	for( size_t y = 0; y < Global::A2powerK[1]; y++ )
        normFactor += n_[0][y];

    float sum = 0.0f;
    // calculate probabilities for order k = 0
    for( size_t y = 0; y < Global::A2powerK[1]; y++ ){
        v_[0][y] = ( ( float )n_[0][y] + A_[0] * 0.25f ) / ( ( float )normFactor + A_[0] );
        sum += v_[0][y];
        range_bar_[0][y] = sum;
    }

    for( size_t k = 1; k < sOrder_+1; k++ ){
        sum = 0.f;
		for( size_t y = 0; y < Global::A2powerK[k+1]; y++ ){
			size_t yk = y / Global::A2powerK[1];
            size_t y2 = y % Global::A2powerK[k];
            v_[k][y] = ( ( float )n_[k][y] + A_[k] * v_[k-1][y2] )/ ( ( float )n_[k-1][yk] + A_[k] );

            if( y % Global::A2powerK[1] == 0 ) sum = 0.f;
            sum += v_[k][y];
            range_bar_[k][y] = sum;
		}
	}
    kmer_freq_is_calculated_ = true;
}

void SeqGenerator::rescale_kmer_frequency( Sequence* refSeq ) {

    size_t L = refSeq->getL();
    // count k-mers for each sequence
    // reset counts for k-mers
    for( size_t k = 0; k < sOrder_+1; k++ ){
        for( size_t y = 0; y < Global::A2powerK[k+1]; y++ ){
            n_seq_[k][y] = 0;
        }
    }
    // count k-mers
    size_t* kmer = refSeq->getKmer();
    for( size_t k = 0; k < sOrder_+1; k++ ){
        // loop over sequence positions
        for( size_t j = k; j < L; j++ ){
            // extract (k+1)-mer
            size_t y = kmer[j] % Global::A2powerK[k+1];
            // count (k+1)mer
            n_seq_[k][y]++;
        }
    }

    // calculate sequence-specific conditional probabilities
    // when k = 0:
    size_t k =0;
    float sum = 0.f;
    for( size_t y = 0; y < Global::A2powerK[k+1]; y++ ){
        v_seq_[k][y] = v_[k][y];
        sum += v_seq_[k][y];
        range_bar_[k][y] = sum;
    }

    // when k = 1:
    k = 1;

    for( size_t y = 0; y < Global::A2powerK[k+1]; y++ ){
        size_t y2 = y % Global::A2powerK[k];
        v_seq_[k][y] = v_[k][y] * ( n_seq_[k][y] + A_[k-1] * v_[k-1][y2] ) / v_[k-1][y2] / ( L + A_[k-1] );
    }
    // re-scale 1st-order background model
    std::vector<float> normFactors( Global::A2powerK[k], 0.0f );
    for( size_t y = 0; y < Global::A2powerK[k+1]; y++ ){
        size_t yk = y / Global::A2powerK[1];
        v_seq_[k][y] = ( n_seq_[k][y] + A_[k] * v_seq_[k][y] ) / ( n_seq_[k-1][yk] + A_[k] );
        normFactors[yk] += v_seq_[k][y];
    }
    // normalization:
    for( size_t y = 0; y < Global::A2powerK[k+1]; y++ ){
        size_t yk = y / Global::A2powerK[1];
        v_seq_[k][y] /= normFactors[yk];
    }

    for( size_t y = 0; y < Global::A2powerK[k+1]; y++ ){
        if( y % Global::A2powerK[1] == 0 )    sum = 0.0f;
        sum += v_seq_[k][y];
        range_bar_[k][y] = sum;
    }

    // for k= 2:
    k = 2;
    // re-scale 2nd-order background model
    for( size_t y = 0; y < Global::A2powerK[k+1]; y++ ){
        size_t y2 = y % Global::A2powerK[k];
        size_t yk = y / Global::A2powerK[1];
        v_seq_[k][y] = ( n_seq_[k][y] + A_[k] * v_seq_[k-1][y2] ) / ( n_seq_[k-1][yk] + A_[k] );
        if( y % Global::A2powerK[1] == 0 ) sum = 0.0f;
        sum += v_seq_[k][y];
        range_bar_[k][y] = sum;
    }

    kmer_freq_is_rescaled_ = true;
}

// generate negative sequences based on k-mer frequencies, given m-fold
std::vector<std::unique_ptr<Sequence>> SeqGenerator::sample_bgseqset_by_fold(size_t fold){

	std::vector<std::unique_ptr<Sequence>> negset;

	calculate_kmer_frequency();
    // todo: can be parallelized
	for( size_t i = 0; i < seqs_.size(); i++ ){
		for( size_t n = 0; n < fold; n++ ){
            if( genericNeg_ ){
                negset.push_back( bg_sequence( seqs_[i]->getL() ) );
            } else {
                negset.push_back( bgseq_on_rescaled_v(seqs_[i]) );
            }
		}
	}

	return negset;

}

// generate negative sequences based on k-mer frequencies, given negative sequence number and max length of input seqs
std::vector<std::unique_ptr<Sequence>> SeqGenerator::sample_bgseqset_by_num(size_t negN, size_t maxL){

    std::vector<std::unique_ptr<Sequence>> negset;

    calculate_kmer_frequency();

    // todo: can be parallelized
    for( size_t n = 0; n < negN; n++ ){
        negset.push_back( bg_sequence(maxL) );
    }

    return negset;
}


// generate each background sequence based on k-mer frequencies from positive set
std::unique_ptr<Sequence> SeqGenerator::bg_sequence(size_t L){

    assert( kmer_freq_is_calculated_ );

    uint8_t* sequence = ( uint8_t* )calloc( L, sizeof( uint8_t ) );
	std::string header = "> bg_seq";

	// sample the first nucleotide
	float random = ( float )rand() / ( float )RAND_MAX;
	for( uint8_t y = 0; y < Global::A2powerK[1]; y++ ){
		if( random <= range_bar_[0][y] ){
			sequence[0] = y+1;
			break;
		}
	}

	// sample the next ( sOrder-1 ) nucleotides
	for( size_t i = 1; i < sOrder_; i++ ){

		// extract y from k-mer
		size_t yk = 0;
		for( size_t k = i; k > 0; k-- ){
			yk += ( sequence[i-k] - 1 ) * Global::A2powerK[k];
		}

		// sample a nucleotide based on k-mer frequency
		random = ( float )rand() / ( float )RAND_MAX;
        for( size_t y = yk, a = 1; y < yk+Global::A2powerK[1]; y++, a++ ){
            sequence[i] = a;
            if( random <= range_bar_[i][y] ){
                break;
            }
        }
	}

    // sample the rest residues from sOrder to L
	for( size_t i = sOrder_; i < L; i++ ){

		// calculate y of K-mer
		size_t yk = 0;
		for( size_t k = sOrder_; k > 0; k-- ){
			yk += ( sequence[i-k] - 1 ) * Global::A2powerK[k];
		}

        random = ( float )rand() / ( float )RAND_MAX;
        for( size_t y = yk, a = 1; y < yk+Global::A2powerK[1]; y++, a++ ){
            sequence[i] = a;
            if( random <= range_bar_[sOrder_][y] ){
			    break;
            }
        }
	}

    std::unique_ptr<Sequence> bg_seq = util::make_unique<Sequence>( sequence, L, header, true );
	if( sequence ) free( sequence );

	return bg_seq;

}

std::unique_ptr<Sequence> SeqGenerator::bgseq_on_rescaled_v( Sequence* refSeq ) {

    assert( kmer_freq_is_calculated_ );

    rescale_kmer_frequency( refSeq );

    assert( kmer_freq_is_rescaled_ );

    size_t L = refSeq->getL();
    uint8_t* sequence = ( uint8_t* )calloc( L, sizeof( uint8_t ) );
    std::string header = "> bg_seq";

    // sample the first nucleotide
    float random = ( float )rand() / ( float )RAND_MAX;
    for( uint8_t y = 0; y < Global::A2powerK[1]; y++ ){
        if( random <= range_bar_[0][y] ){
            sequence[0] = y+1;
            break;
        }
    }

    // sample the next ( sOrder-1 ) nucleotides
    for( size_t i = 1; i < sOrder_; i++ ){

        // extract y from k-mer
        size_t yk = 0;
        for( size_t k = i; k > 0; k-- ){
            yk += ( sequence[i-k] - 1 ) * Global::A2powerK[k];
        }

        // sample a nucleotide based on k-mer frequency
        random = ( float )rand() / ( float )RAND_MAX;
        for( size_t y = yk, a = 1; y < yk+Global::A2powerK[1]; y++, a++ ){
            sequence[i] = a;
            if( random <= range_bar_[i][y] ){
                break;
            }
        }
    }

    // sample the rest residues from sOrder to L
    for( size_t i = sOrder_; i < L; i++ ){

        // calculate y of K-mer
        size_t yk = 0;
        for( size_t k = sOrder_; k > 0; k-- ){
            yk += ( sequence[i-k] - 1 ) * Global::A2powerK[k];
        }

        random = ( float )rand() / ( float )RAND_MAX;
        for( size_t y = yk, a = 1; y < yk+Global::A2powerK[1]; y++, a++ ){
            sequence[i] = a;
            if( random <= range_bar_[sOrder_][y] ){
                break;
            }
        }
    }

    std::unique_ptr<Sequence> bg_seq = util::make_unique<Sequence>( sequence, L, header, true );
    if( sequence ) free( sequence );

    return bg_seq;

}

// copy sequences from positive set
std::unique_ptr<Sequence> SeqGenerator::raw_sequence( Sequence* seq ){

    size_t L = seq->getL();
    uint8_t* sequence = ( uint8_t* )calloc( L, sizeof( uint8_t ) );
    std::string header = "> raw_seq";

    // copy original sequence to the raw sequence
    for( size_t i = 0; i < L; i++ ){
        sequence[i] = seq->getSequence()[i];
    }

    std::unique_ptr<Sequence> raw_seq = util::make_unique<Sequence>( sequence, L, header, true );

    if( sequence ) free( sequence );

    return raw_seq;

}

// generate sequences with given motif embedded into the given sequences
std::vector<std::unique_ptr<Sequence>> SeqGenerator::arti_posset_motif_embedded( size_t at ){

	std::vector<std::unique_ptr<Sequence>> posset_with_motif_embedded;

	calculate_kmer_frequency();

    size_t seq_size = (size_t)( (float)seqs_.size() * q_ );

    // generative sequence with motif embedded for the q portion
	for( size_t i = 0; i < seq_size; i++ ){
        posset_with_motif_embedded.push_back(posseq_motif_embedded( seqs_[i], at ));
	}
    // copy the original sequences for the rest 1-q portion
    for( size_t i = seq_size; i < seqs_.size(); i++ ){
        posset_with_motif_embedded.push_back( raw_sequence( seqs_[i] ) );
    }

    // randomly shuffle the sequence set after implantation
    std::random_shuffle( posset_with_motif_embedded.begin(), posset_with_motif_embedded.end() );

    return posset_with_motif_embedded;
}

// generate sequences with given motif embedded into the given sequences
std::unique_ptr<Sequence> SeqGenerator::posseq_motif_embedded( Sequence* seq, size_t at ){

    size_t W = motif_->getW();
    size_t L = seq -> getL();
    size_t LW = L + W;
    uint8_t* sequence = ( uint8_t* )calloc( LW, sizeof( uint8_t ) );
    std::string header = seq->getHeader() + "+motif";
    // sample nucleotides till the position for motif implantation
    // due to k-mer frequencies
    if( at == 0 ) {
        // implant motif due to uniform distribution
        std::uniform_int_distribution<> range((int)sOrder_, (int)(L-W+1));
        at = (size_t) range(rngx_);
    } else{
        // implant motif around the given position due to normal distribution
        // set mean as the given position and variance as 10
        std::normal_distribution<> nd{(float)at, 10.f};
        at = (size_t) std::round( nd(rngx_) );
        at = (at > 0) ? at : 0;
        at = (at < L) ? at : L;
    }

    // copy the left part of the given sequence
    for( size_t i = 0; i < at; i++ ){
        sequence[i] = seq->getSequence()[i];
    }

    // insert motif in the middle
    // sample W-nt based on conditional probabilities of the motif
    for( size_t j = 0; j < W; j++ ){
        // extract y from K-mer
        size_t yk = 0;
        for( size_t k = sOrder_; k > 0; k-- ){
            yk += ( sequence[j+at-k] - 1 ) * Global::A2powerK[k];
        }
        // sample a nucleotide based on k-mer frequency
        float random = ( float )rand() / ( float )RAND_MAX;
        float f = 0.0f;
        for( uint8_t a = 1; a <= Global::A2powerK[1]; a++ ){
            f += motif_->getV()[sOrder_][yk+a-1][j];
            sequence[j+at] = a;
            if( random <= f ){
                break;
            }
        }
    }

    // copy the right part of the given sequence
    for( size_t i = at+W; i < L; i++ ){
        sequence[i] = seq->getSequence()[i-W];
    }

    std::unique_ptr<Sequence> seq_with_motif = util::make_unique<Sequence>( sequence, L, header, true );

    return seq_with_motif;

}

std::vector<std::unique_ptr<Sequence>> SeqGenerator::seqset_with_motif_masked(float **r){

	std::vector<std::unique_ptr<Sequence>> seqset;
	size_t W = motif_->getW();

	for( size_t n = 0; n < seqs_.size(); n++ ){
		seqset.push_back(sequence_with_motif_masked( seqs_[n], W, r[n]) );
	}

	return seqset;
}

std::unique_ptr<Sequence> SeqGenerator::sequence_with_motif_masked(Sequence *posseq, size_t W, float *r){

	float cutoff = 0.3f;
    size_t L = posseq->getL();

    uint8_t* masked_seq = ( uint8_t* )calloc( L, sizeof( uint8_t ) );
	// copy the original sequence from a C-array to a vector
	// and mask patterns with a weight larger than certain cutoff
	size_t i = 0;
    size_t j = 0;
	while( i <= L - W ){
		if( r[L-W-i] < cutoff ){
            masked_seq[j] = posseq->getSequence()[i];
			i++;
            j++;
		} else {
			i += W;
		}
	}

	std::string header = "> seq with motif masked";

    std::unique_ptr<Sequence> seq_mask_motif = util::make_unique<Sequence>( masked_seq, L, header, true );

	return seq_mask_motif;
}

void SeqGenerator::write( char* odir, std::string basename, std::vector<std::unique_ptr<Sequence>> seqset ){
	/**
	 * save the generated sequence set in fasta file:
	 */
	std::string opath = std::string( odir ) + '/' + basename + ".fasta";

	std::ofstream ofile( opath );

	for( size_t n = 0; n < seqset.size(); n++ ){
		ofile << seqset[n]->getHeader() << std::endl;
		for( size_t i = 0; i < seqset[n]->getL(); i++ ){
			ofile << Alphabet::getBase( seqset[n]->getSequence()[i] );
		}
		ofile << std::endl;
	}

}

void SeqGenerator::write( char* odir, std::string basename, std::vector<Sequence*> seqset ){
    /**
     * save the generated sequence set in fasta file:
     */
    std::string opath = std::string( odir ) + '/' + basename + ".fasta";

    std::ofstream ofile( opath );

    for( size_t n = 0; n < seqset.size(); n++ ){
        ofile << seqset[n]->getHeader() << std::endl;
        for( size_t i = 0; i < seqset[n]->getL(); i++ ){
            ofile << Alphabet::getBase( seqset[n]->getSequence()[i] );
        }
        ofile << std::endl;
    }

}

