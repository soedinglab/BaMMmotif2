#include "SeqGenerator.h"

SeqGenerator::SeqGenerator( std::vector<Sequence*> seqs, Motif* motif, size_t sOrder, float q ){

	seqs_ = seqs;
	sOrder_ = sOrder;
    motif_ = motif;
    q_ = q;

	for( size_t k = 0; k < sOrder_ + 8; k++ ){
		Y_.push_back( ipow( Alphabet::getSize(), k ) );
	}

	freqs_ = ( float** )calloc( sOrder_+1, sizeof( float* ) );
	count_ = ( size_t** )calloc( sOrder_+1, sizeof( size_t* ) );
	for( size_t k = 0; k < sOrder_+1; k++ ){
		freqs_[k] = ( float* )calloc( Y_[k+1], sizeof( float ) );
		count_[k] = ( size_t* )calloc( Y_[k+1], sizeof( size_t ) );
	}

}

SeqGenerator::~SeqGenerator(){

	for( size_t k = 0; k < sOrder_+1; k++ ){
		free( freqs_[k] );
		free( count_[k] );
	}
	free( freqs_ );
	free( count_ );

}

void SeqGenerator::calculate_kmer_frequency(){

	// reset counts for k-mers
	for( size_t k = 0; k < sOrder_+1; k++ ){
		for( size_t y = 0; y < Y_[k+1]; y++ ){
			count_[k][y] = 0;
		}
	}

	// count k-mers
	for( size_t i = 0; i < seqs_.size(); i++ ){
		size_t L = seqs_[i]->getL();
		size_t* kmer = seqs_[i]->getKmer();
		for( size_t k = 0; k < sOrder_+1; k++ ){
			// loop over sequence positions
			for( size_t j = k; j < L; j++ ){
				// extract (k+1)-mer
				size_t y = kmer[j] % Y_[k+1];
				// count (k+1)mer
				count_[k][y]++;
			}
		}
	}

	// calculate frequencies
	size_t normFactor = 0;
	for( size_t y = 0; y < Y_[1]; y++ )	normFactor += count_[0][y];
	for( size_t y = 0; y < Y_[1]; y++ )	freqs_[0][y] = ( float )count_[0][y] / ( float )normFactor;
	for( size_t k = 1; k < sOrder_+1; k++ ){
		for( size_t y = 0; y < Y_[k+1]; y++ ){
			size_t yk = y / Y_[1];
			freqs_[k][y] = ( float )count_[k][y] / ( float )count_[k-1][yk];
		}
	}
}


// generate negative sequences based on k-mer frequencies from positive set
std::vector<std::unique_ptr<Sequence>> SeqGenerator::arti_negset(
		size_t fold ){

	std::vector<std::unique_ptr<Sequence>> negset;

	calculate_kmer_frequency();

	for( size_t i = 0; i < seqs_.size(); i++ ){
		size_t L = seqs_[i]->getL();
		for( size_t n = 0; n < fold; n++ ){
			negset.push_back( negseq_dimer_freq( L ) );
		}
	}
	return negset;
}

// generate background sequence based on k-mer frequencies from positive set
std::unique_ptr<Sequence> SeqGenerator::negseq_dimer_freq( size_t L ){

	uint8_t* sequence = ( uint8_t* )calloc( L, sizeof( uint8_t ) );
	std::string header = "negative sequence";

	// sample the first nucleotide
	double random = ( double )rand() / ( double )RAND_MAX;
	double f = 0.0f;
	for( uint8_t a = 1; a <= Y_[1]; a++ ){
		f += freqs_[0][a-1];
		if( random <= f ){
			sequence[0] = a;
			break;
		}
		// Trick: this solves the numerical problem
		if( sequence[0] == 0 )	sequence[0] = a;
	}

	// sample the next ( sOrder-1 ) nucleotides
	for( size_t i = 1; i < sOrder_; i++ ){

		// extract y from k-mer
		size_t yk = 0;
		for( size_t k = i; k > 0; k-- ){
			yk += ( sequence[i-k] - 1 ) * Y_[k];
		}

		// sample a nucleotide based on k-mer frequency
		random = ( double )rand() / ( double )RAND_MAX;
		f = 0.0f;
		for( uint8_t a = 1; a <= Y_[1]; a++ ){
			f += freqs_[i][yk+a-1];
			if( random <= f ){
				sequence[i] = a;
				break;
			}
			// Trick: this solves the numerical problem
			if( sequence[i] == 0 )	sequence[i] = a;
		}
	}

	for( size_t i = sOrder_; i < L; i++ ){
		random = ( double )rand() / ( double )RAND_MAX;
		// calculate y of K-mer
		size_t yk = 0;
		for( size_t k = sOrder_; k > 0; k-- ){
			yk += ( sequence[i-k] - 1 ) * Y_[k];
		}

		// assign a nucleotide based on k-mer frequency
		f = 0.0f;
		for( uint8_t a = 1; a <= Y_[1]; a++ ){
			f += freqs_[sOrder_][yk+a-1];
			if( random <= f ){
				sequence[i] = a;
				break;
			}
			// Trick: this solves the numerical problem
			if( sequence[i] == 0 )	sequence[i] = a;
		}
	}

	std::unique_ptr<Sequence> seq( new Sequence( sequence, L, header, Y_ ) );

	free( sequence );

	return seq;

}

// generate pseudo-positive sequences based on each test set
std::vector<std::unique_ptr<Sequence>> SeqGenerator::arti_posset_motif_embedded( size_t fold ){

	std::vector<std::unique_ptr<Sequence>> posset_with_motif_embedded;

	calculate_kmer_frequency();

	for( size_t i = 0; i < seqs_.size(); i++ ){
		size_t L = seqs_[i]->getL();
		for( size_t n = 0; n < fold; n++ ){
			posset_with_motif_embedded.push_back(posseq_motif_embedded(L) );
		}
	}
	return posset_with_motif_embedded;
}

// generate pseudo-sequence based on k-mer frequencies from positive set
std::unique_ptr<Sequence> SeqGenerator::posseq_motif_embedded(size_t L){

	uint8_t* sequence = ( uint8_t* )calloc( L, sizeof( uint8_t ) );

	std::string header = "artificial sequence with motif embedded in the middle";

	for( uint8_t i = 0; i < L; i++ ){
		sequence[i] = 1;
	}

	// sample the first nucleotide
	double random = ( double )rand() / ( double )RAND_MAX;
	double f = 0.0f;
	for( uint8_t a = 1; a <= Y_[1]; a++ ){
		f += freqs_[0][a-1];
		if( random <= f ){
			sequence[0] = a;
			break;
		}
		// Trick: this solves the numerical problem
		if( sequence[0] == 0 )	sequence[0] = a;
	}

	// sample the next ( K-1 ) nucleotides
	for( size_t i = 1; i < sOrder_; i++ ){

		// extract y from k-mer
		size_t yk = 0;
		for( size_t k = i; k > 0; k-- ){
			yk += ( sequence[i-k] - 1 ) * Y_[k];
		}

		// sample a nucleotide based on k-mer frequency
		random = ( double )rand() / ( double )RAND_MAX;
		f = 0.0f;
		for( uint8_t a = 1; a <= Y_[1]; a++ ){
			f += freqs_[i][yk+a-1];
			if( random <= f ){
				sequence[i] = a;
				break;
			}
			// Trick: this solves the numerical problem
			if( sequence[i] == 0 )	sequence[i] = a;
		}
	}
	// sample nucleotides till the half length of the sequence
	// due to k-mer frequencies
	size_t mid = L / 2;
	for( size_t i = sOrder_; i < mid; i++ ){

		// extract y from K-mer
		size_t yk = 0;
		for( size_t k = sOrder_; k > 0; k-- ){
			yk += ( sequence[i-k] - 1 ) * Y_[k];
		}

		// sample a nucleotide based on k-mer frequency
		random = ( double )rand() / ( double )RAND_MAX;
		f = 0.0f;
		for( uint8_t a = 1; a <= Y_[1]; a++ ){
			f += freqs_[sOrder_][yk+a-1];
			if( random <= f ){
				sequence[i] = a;
				break;
			}
			// Trick: this solves the numerical problem
			if( sequence[i] == 0 )	sequence[i] = a;
		}
	}

	// sample W-nt based on conditional probabilities of the motif
	size_t W = motif_->getW();
	float*** v = motif_->getV();
	for( size_t j = 0; j < W; j++ ){
		// extract y from K-mer
		size_t yk = 0;
		for( size_t k = sOrder_; k > 0; k-- ){
			yk += ( sequence[j+mid-k] - 1 ) * Y_[k];
		}
		// sample a nucleotide based on k-mer frequency
		random = ( double )rand() / ( double )RAND_MAX;
		f = 0.0f;
		for( uint8_t a = 1; a <= Y_[1]; a++ ){
			f += v[sOrder_][yk+a-1][j];
			if( random <= f ){
				sequence[j+mid] = a;
				break;
			}
			// Trick: this solves the numerical problem
			if( sequence[j+mid] == 0 )	sequence[j+mid] = a;
		}

	}

	// sample the rest of the sequence till it reaches the length of L
	for( size_t i = mid + W; i < L; i++ ){
		// extract y from K-mer
		size_t yk = 0;
		for( size_t k = sOrder_; k > 0; k-- ){
			yk += ( sequence[i-k] - 1 ) * Y_[k];
		}

		// sample a nucleotide based on k-mer frequency
		random = ( double )rand() / ( double )RAND_MAX;
		f = 0.0f;
		for( uint8_t a = 1; a <= Y_[1]; a++ ){
			f += freqs_[sOrder_][yk+a-1];
			if( random <= f ){
				sequence[i] = a;
				break;
			}
			// Trick: this solves the numerical problem
			if( sequence[i] == 0 )	sequence[i] = a;
		}
	}

	std::unique_ptr<Sequence> seq( new Sequence( sequence, L, header, Y_ ) );

	free( sequence );

	return seq;

}

std::vector<std::unique_ptr<Sequence>> SeqGenerator::arti_negset_motif_masked( float** r ){

	std::vector<std::unique_ptr<Sequence>> seqset;
	size_t W = motif_->getW();

	for( size_t n = 0; n < seqs_.size(); n++ ){

		seqset.push_back( negseq_motif_masked( seqs_[n], W, r[n] ) );
	}

	return seqset;
}

std::unique_ptr<Sequence> SeqGenerator::negseq_motif_masked( Sequence* posseq, size_t W, float* r  ){

	float cutoff = 0.3f;

	std::vector<uint8_t> fake_seq;
	// copy the original sequence from a C-array to a vector
	// and mask patterns with a weight larger than certain cutoff
	size_t i = 0;
	while( i < posseq->getL() ){
		size_t LW1 = posseq->getL() - W + 1;
		if( r[LW1-i+1] < cutoff ){
			fake_seq.push_back( posseq->getSequence()[i] );
			i++;
		} else {
			i += W;
		}
	}

	size_t L = fake_seq.size();
	uint8_t* real_seq = ( uint8_t* )calloc( L, sizeof( uint8_t ) );
	for( size_t j = 0; j < L; j++ ){
		real_seq[j] = fake_seq[j];
	}
	std::string header = "sequence with masked motif";
	std::unique_ptr<Sequence> seq_mask_motif( new Sequence( real_seq, L, header, Y_, true ) );
	free( real_seq );
	return seq_mask_motif;
}

void SeqGenerator::write( char* odir, std::string basename, std::vector<std::unique_ptr<Sequence>> seqset ){
	/**
	 * save the generated sequence set in fasta file:
	 * (1) posSequenceBasename_pseudo.fasta
	 */

	std::string opath = std::string( odir ) + '/' + basename + ".fasta";

	std::ofstream ofile( opath );

	for( size_t n = 0; n < seqset.size(); n++ ){
		ofile << "> " << seqset[n]->getHeader() << std::endl;
		for( size_t i = 0; i < seqset[n]->getL(); i++ ){
			ofile << Alphabet::getBase( seqset[n]->getSequence()[i] );
		}
		ofile << std::endl;
	}

}

