#include "SeqGenerator.h"

SeqGenerator::SeqGenerator( std::vector<Sequence*> seqs, Motif* motif, int sOrder ){

	seqs_ = seqs;
	motif_= motif;
	sOrder_ = sOrder;

	for( int k = 0; k < sOrder_ + 4; k++ ){
		Y_.push_back( ipow( Alphabet::getSize(), k ) );
	}

	freqs_ = ( float** )calloc( sOrder_+1, sizeof( float* ) );
	count_ = ( int** )calloc( sOrder_+1, sizeof( int* ) );
	for( int k = 0; k < sOrder_+1; k++ ){
		freqs_[k] = ( float* )calloc( Y_[k+1], sizeof( float ) );
		count_[k] = ( int* )calloc( Y_[k+1], sizeof( int ) );
	}

}

SeqGenerator::~SeqGenerator(){

	for( int k = 0; k < sOrder_+1; k++ ){
		free( freqs_[k] );
		free( count_[k] );
	}
	free( freqs_ );
	free( count_ );

}

void SeqGenerator::calculate_kmer_frequency(){

	// reset counts for k-mers
	for( int k = 0; k < sOrder_+1; k++ ){
		for( int y = 0; y < Y_[k+1]; y++ ){
			count_[k][y] = 0;
		}
	}

	// count k-mers
	for( size_t i = 0; i < seqs_.size(); i++ ){
		int L = seqs_[i]->getL();
		int* kmer = seqs_[i]->getKmer();
		for( int k = 0; k < sOrder_+1; k++ ){
			// loop over sequence positions
			for( int j = k; j < L; j++ ){
				// extract (k+1)-mer
				int y = kmer[j] % Y_[k+1];
				// skip non-defined alphabet letters
				if( y >= 0 ){
					// count (k+1)mer
					count_[k][y]++;
				}
			}
		}
	}

	// calculate frequencies
	int normFactor = 0;
	for( int y = 0; y < Y_[1]; y++ )	normFactor += count_[0][y];
	for( int y = 0; y < Y_[1]; y++ )	freqs_[0][y] = ( float )count_[0][y] / ( float )normFactor;
	for( int k = 1; k < sOrder_+1; k++ ){
		for( int y = 0; y < Y_[k+1]; y++ ){
			int yk = y / Y_[1];
			freqs_[k][y] = ( float )count_[k][y] / ( float )count_[k-1][yk];
		}
	}

}


// generate negative sequences
std::vector<std::unique_ptr<Sequence>> SeqGenerator::sample_negative_seqset( int fold ){

	std::vector<std::unique_ptr<Sequence>> sampleSet;

	calculate_kmer_frequency();

	for( size_t i = 0; i < seqs_.size(); i++ ){
		int L = seqs_[i]->getL();
		for( int n = 0; n < fold; n++ ){
			sampleSet.push_back( sample_negative_sequence( L ) );
		}
	}
	return sampleSet;
}

// generate background sequence based on k-mer frequencies from positive set
std::unique_ptr<Sequence> SeqGenerator::sample_negative_sequence( int L ){

	uint8_t* sequence = ( uint8_t* )calloc( L, sizeof( uint8_t ) );
	std::string header = "> negative sequence";

	// sample the first nucleotide
	double random = ( double )rand() / ( double )RAND_MAX;
	double f = 0.0f;
	for( uint8_t a = 1; a <= Y_[1]; a++ ){
		f += freqs_[0][a-1];
		if( random <= f ){
			sequence[0] = a;
			break;
		}
		if( sequence[0] == 0 )	sequence[0] = a;		// Trick: this is to solve the numerical problem
	}

	// sample the next ( sOrder-1 ) nucleotides
	for( int i = 1; i < sOrder_; i++ ){

		// extract y from k-mer
		int yk = 0;
		for( int k = i; k > 0; k-- ){
			yk += ( sequence[i-k] - 1 ) * Y_[k];
		}

		// sample a nucleotide based on k-mer frequency
		random = ( double )rand() / ( double )RAND_MAX;	// get another random double number
		f = 0.0f;
		for( uint8_t a = 1; a <= Y_[1]; a++ ){
			f += freqs_[i][yk+a-1];
			if( random <= f ){
				sequence[i] = a;
				break;
			}
			if( sequence[i] == 0 )	sequence[i] = a;	// Trick: this is to solve the numerical problem
		}
	}

	for( int i = sOrder_; i < L; i++ ){
		random = ( double )rand() / ( double )RAND_MAX;	// get another random double number
		// calculate y of K-mer
		int yk = 0;
		for( int k = sOrder_; k > 0; k-- ){
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
			if( sequence[i] == 0 )	sequence[i] = a;	// Trick: this is to solve the numerical problem
		}
	}

	std::unique_ptr<Sequence> seq( new Sequence( sequence, L, header, Y_ ) );

	free( sequence );

	return seq;

}

// generate pseudo-positive sequences based on each test set
std::vector<std::unique_ptr<Sequence>> SeqGenerator::sample_pseudo_seqset( int fold ){

	std::vector<std::unique_ptr<Sequence>> sampleSet;

	calculate_kmer_frequency();

	for( size_t i = 0; i < seqs_.size(); i++ ){
		int L = seqs_[i]->getL();
		for( int n = 0; n < fold; n++ ){
			sampleSet.push_back( sample_pseudo_sequence( L ) );
		}
	}
	return sampleSet;
}

// generate pseudo-sequence based on k-mer frequencies from positive set
std::unique_ptr<Sequence> SeqGenerator::sample_pseudo_sequence( int L ){

	uint8_t a;
	int i,j, k, yk;

	uint8_t* sequence = ( uint8_t* )calloc( L, sizeof( uint8_t ) );
	std::string header = "> pseudo sequence";

	for( i = 0; i < L; i++ ){
		sequence[i] = 1;
	}

	// sample the first nucleotide
	double random = ( double )rand() / ( double )RAND_MAX;
	double f = 0.0f;
	for( a = 1; a <= Y_[1]; a++ ){
		f += freqs_[0][a-1];
		if( random <= f ){
			sequence[0] = a;
			break;
		}
		if( sequence[0] == 0 )	sequence[0] = a;		// Trick: this is to solve the numerical problem
	}

	// sample the next ( K-1 ) nucleotides
	for( i = 1; i < sOrder_; i++ ){

		// extract y from k-mer
		yk = 0;
		for( k = i; k > 0; k-- ){
			yk += ( sequence[i-k] - 1 ) * Y_[k];
		}

		// sample a nucleotide based on k-mer frequency
		random = ( double )rand() / ( double )RAND_MAX;	// get another random double number
		f = 0.0f;
		for( a = 1; a <= Y_[1]; a++ ){
			f += freqs_[i][yk+a-1];
			if( random <= f ){
				sequence[i] = a;
				break;
			}
			if( sequence[i] == 0 )	sequence[i] = a;	// Trick: this is to solve the numerical problem
		}
	}
	// sample nucleotides till the half length of the sequence, due to k-mer frequencies
	int mid = L / 2;
	for( i = sOrder_; i < mid; i++ ){

		// extract y from K-mer
		yk = 0;
		for( k = sOrder_; k > 0; k-- ){
			yk += ( sequence[i-k] - 1 ) * Y_[k];
		}

		// sample a nucleotide based on k-mer frequency
		random = ( double )rand() / ( double )RAND_MAX;	// get another random double number
		f = 0.0f;
		for( a = 1; a <= Y_[1]; a++ ){
			f += freqs_[sOrder_][yk+a-1];
			if( random <= f ){
				sequence[i] = a;
				break;
			}
			if( sequence[i] == 0 )	sequence[i] = a;	// Trick: this is to solve the numerical problem
		}
	}

	// sample W-length nucleotides based on conditional probabilities of the motif
	int W = motif_->getW();
	float*** v = motif_->getV();
	for( j = 0; j < W; j++ ){
		// extract y from K-mer
		yk = 0;
		for( k = sOrder_; k > 0; k-- ){
			yk += ( sequence[j+mid-k] - 1 ) * Y_[k];
		}
		// sample a nucleotide based on k-mer frequency
		random = ( double )rand() / ( double )RAND_MAX;	// get another random double number
		f = 0.0f;
		for( a = 1; a <= Y_[1]; a++ ){
			f += v[sOrder_][yk+a-1][j];
			if( random <= f ){
				sequence[j+mid] = a;
				break;
			}
			if( sequence[j+mid] == 0 )	sequence[j+mid] = a;	// Trick: this is to solve the numerical problem
		}

	}

	// sample the rest of the sequence till it reaches the length of L
	for( i = mid + W; i < L; i++ ){
		// extract y from K-mer
		yk = 0;
		for( k = sOrder_; k > 0; k-- ){
			yk += ( sequence[i-k] - 1 ) * Y_[k];
		}

		// sample a nucleotide based on k-mer frequency
		random = ( double )rand() / ( double )RAND_MAX;	// get another random double number
		f = 0.0f;
		for( a = 1; a <= Y_[1]; a++ ){
			f += freqs_[sOrder_][yk+a-1];
			if( random <= f ){
				sequence[i] = a;
				break;
			}
			if( sequence[i] == 0 )	sequence[i] = a;	// Trick: this is to solve the numerical problem
		}
	}

	std::unique_ptr<Sequence> seq( new Sequence( sequence, L, header, Y_ ) );

	free( sequence );

	return seq;

}

void SeqGenerator::write_pseudoset(){

	/**
	 * save the generated sequence set in one file:
	 * (1) posSequenceBasename_pseudo.fasta
	 */

	std::string opath = std::string( Global::outputDirectory ) + '/'
			+ Global::posSequenceBasename + "_pseudo.fasta";

	std::ofstream ofile( opath );

	std::vector<std::unique_ptr<Sequence>> seqs = sample_pseudo_seqset( Global::mFold );
	for( size_t n = 0; n < seqs.size(); n++ ){
		ofile << "> " << seqs[n]->getHeader() << std::endl;
		for( int i = 0; i < seqs[n]->getL(); i++ ){
			ofile << Alphabet::getBase( seqs[n]->getSequence()[i] );
		}
		ofile << std::endl;
	}

}
