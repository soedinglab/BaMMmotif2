/*
 * BackgroundModel.cpp
 *
 *  Created on: Apr 18, 2016
 *      Author: administrator
 */

#include "BackgroundModel.h"
#include "SequenceSet.h"

BackgroundModel::BackgroundModel(){
	powA_ = new int[Global::bgModelOrder+2];
	for( int k = 0; k < Global::bgModelOrder + 2; k++ )
		powA_[k] = Global::ipow( Alphabet::getSize(), k );

	Y_ = 0;
	for( int k = 0; k < Global::bgModelOrder + 1; k++ )
		Y_ += powA_[k+1];

	n_occ_ = ( int* )calloc( Y_, sizeof( int ) );

	v_ = ( float* )calloc( Y_, sizeof( float ) );
}

BackgroundModel::~BackgroundModel(){
	free( v_ );
	free( n_occ_ );
	delete[] powA_;
}

void BackgroundModel::init( std::vector<int> folds = std::vector<int> () ){

	// count n_occ[y]
	if( folds.size() == 0 ){
		// when the full negSeqSet is applied
		for( unsigned int n = 0; n < Global::negSequenceSet->getN(); n++ ){
			Sequence seq = Global::negSequenceSet->getSequences()[n];
			for( int k = 0; k < Global::bgModelOrder + 1; k++ ){
				for( int i = k; i < seq.getL(); i++ ){
					int y = seq.extractKmerbg( i, k );
					if( y >= 0 )				// skip the non-set alphabets, such as N
						n_occ_[y]++;
				}
			}
		}
	} else {
		// for training sets in cross-validation
		for( unsigned int f = 0; f < folds.size(); f++ ){
			for( unsigned int r = 0; r < Global::negFoldIndices[f].size(); r++ ){
				unsigned int idx = Global::negFoldIndices[f][r];
				Sequence seq = Global::negSequenceSet->getSequences()[idx];
				for( int k = 0; k < Global::bgModelOrder + 1; k++ ){
					for( int i = k; i < seq.getL(); i++ ){
						int y = seq.extractKmerbg( i, k );
						if( y >= 0 )			// skip the non-set alphabets, such as N
							n_occ_[y]++;
					}
				}
			}
		}
	}

	// Only for testing
//	if( Global::verbose ){
//		printf( " ___________________________\n"
//				"|                           |\n"
//				"| k-mer counts for bgModel  |\n"
//				"|___________________________|\n\n" );
//		for( int y = 0; y < Y_; y++ )
//			std::cout << n_occ_[y] << '\t';
//		std::cout << std::endl;
//	}

    // calculate v from k-mer counts n_occ
	calculateV();

//	if( Global::verbose ) print();
}

void BackgroundModel::calculateV(){

	// sum of the counts n[y]
	float sum = 0.0f;
	for( int y = 0; y < powA_[1]; y++ )
		sum += float( n_occ_[y] );

	// for k = 0: v = freqs
	for( int y = 0; y < powA_[1]; y++ )
		v_[y] = n_occ_[y] / sum;

	// for k > 0:
	for( int y = powA_[1]; y < Y_; y++ ){
		int y2 = y / powA_[1] - 1;						// cut off the first nucleotide, e.g. from (5'-)ACGT(-3') to CGT
		for( int k = 1; k < Global::bgModelOrder; k++ ){
			int yk = y % powA_[k];						// cut off the last nucleotide, e.g. from (5'-)ACGT(-3') to ACG
			v_[y] = ( n_occ_[y] + Global::bgModelAlpha * v_[y2] ) / ( n_occ_[yk] + Global::bgModelAlpha );
//			v_[y] =  float( n_occ[y] ) / n_occ[yk];		// without interpolation of alpha
		}
	}
}

float* BackgroundModel::getV(){
    return v_;
}

void BackgroundModel::print(){
	printf( " ___________________________\n"
			"|                           |\n"
			"| probabilities for bgModel |\n"
			"|___________________________|\n\n" );
	for( int k = 0, y = 0; k < Global::bgModelOrder + 1; k ++ ){
		for( int i = 0; i < powA_[k+1]; i++, y++ )
			std::cout << std::scientific << v_[y] << '\t';
		std::cout << std::endl;
	}
}

void BackgroundModel::write(){
	std::string opath = std::string( Global::outputDirectory )  + '/' + std::string( Global::posSequenceBasename ) + ".probsBg";
	std::ofstream ofile( opath.c_str() );
	for( int k = 0, y = 0; k < Global::bgModelOrder + 1; k++ ){
		for( int i = 0; i < powA_[k+1]; i++, y++ )
			ofile << std::scientific << v_[y] << '\t';
		ofile << std::endl;
	}
}
