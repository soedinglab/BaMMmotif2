/*
 * BackgroundModel.cpp
 *
 *  Created on: Apr 18, 2016
 *      Author: administrator
 */

#include "BackgroundModel.h"
#include "SequenceSet.h"

BackgroundModel::BackgroundModel(){
	n_bg_ = ( int** )calloc( K_+1, sizeof( int* ) );
	v_bg_ = ( float** )calloc( K_+1, sizeof( float* ) );
	for( int k = 0; k < K_+1; k++ ){
		n_bg_[k] = ( int* )calloc( Global::powA[k+1], sizeof( int ) );
		v_bg_[k] = ( float* )calloc( Global::powA[k+1], sizeof( float ) );
	}
}

BackgroundModel::~BackgroundModel(){
	for( int k = 0; k < K_+1; k++ ){
		free( v_bg_[k] );
		free( n_bg_[k] );
	}
	free( v_bg_ );
	free( n_bg_ );
}

void BackgroundModel::init( std::vector<int> folds = std::vector<int> () ){

	// count kmers
	if( folds.size() == 0 ){
		// when the full negSeqSet is applied
		for( unsigned int n = 0; n < Global::negSequenceSet->getN(); n++ ){
			Sequence seq = Global::negSequenceSet->getSequences()[n];
			for( int k = 0; k < K_ + 1; k++ ){
				for( int i = 0; i < seq.getL(); i++ ){
					int y = seq.extractKmer( i, k );
					if( y >= 0 )				// skip the non-set alphabets, such as N
						n_bg_[k][y]++;
				}
			}
		}
	} else {
		// for training sets in cross-validation
		for( unsigned int f = 0; f < folds.size(); f++ ){
			for( unsigned int r = 0; r < Global::negFoldIndices[f].size(); r++ ){
				unsigned int idx = Global::negFoldIndices[f][r];
				Sequence seq = Global::negSequenceSet->getSequences()[idx];
				for( int k = 0; k < K_ + 1; k++ ){
					for( int i = k; i < seq.getL(); i++ ){
						int y = seq.extractKmer( i, k );
						if( y >= 0 )			// skip the non-set alphabets, such as N
							n_bg_[k][y]++;
					}
				}
			}
		}
	}

//	 // Only for testing
//	if( Global::verbose ){
//		printf( " ___________________________\n"
//				"|                           |\n"
//				"| k-mer counts for bgModel  |\n"
//				"|___________________________|\n\n" );
//		for( int k = 0; k < K_+1; k++ ){
//			for( int y = 0; y < powA_[k+1]; y++ )
//				std::cout << n_bg_[k][y] << '\t';
//			std::cout << std::endl;
//		}
//	}

    // calculate v from k-mer counts n_occ
	calculateVbg();

	if( Global::verbose ) print();

	write();
}

void BackgroundModel::calculateVbg(){

	// sum over counts for all mono-nucleotides
	float sum = 0.0f;
	for( int y = 0; y < Global::powA[1]; y++ )
		sum += float( n_bg_[0][y] );

	// for k = 0: v = frequency
	for( int y = 0; y < Global::powA[1]; y++ )
		v_bg_[0][y] = n_bg_[0][y] / sum;

	// for k > 0:
	int y2;										// cut off the first nucleotide in (k+1)-mer y, e.g. from ACGT to CGT
	int yk;										// cut off the last nucleotide in (k+1)-mer y , e.g. from ACGT to ACG
	for( int k = 1; k < K_+1; k++ ){
		for( int y = 0; y < Global::powA[k+1]; y++ ){
			y2 = y % Global::powA[k];
			yk = y / Global::powA[1];
			v_bg_[k][y] = ( n_bg_[k][y] + Global::bgModelAlpha * v_bg_[k-1][y2] )
						/ ( n_bg_[k-1][yk] + Global::bgModelAlpha );
		}
	}
}

float** BackgroundModel::getVbg(){
    return v_bg_;
}

void BackgroundModel::print(){
	printf( " ___________________________\n"
			"|                           |\n"
			"| probabilities for bgModel |\n"
			"|___________________________|\n\n" );
	for( int k = 0; k < K_+1; k++ ){
		for( int y = 0; y < Global::powA[k+1]; y++ )
			std::cout << std::fixed << std::setprecision(6) << v_bg_[k][y] << '\t';
		std::cout << std::endl;
	}
}

void BackgroundModel::write(){
	std::string opath = std::string( Global::outputDirectory )  + '/'
			+ std::string( Global::posSequenceBasename ) + ".condsBg";
	std::ofstream ofile( opath.c_str() );
	for( int k = 0; k < K_+1; k++ ){
		for( int y = 0; y < Global::powA[k+1]; y++ )
			ofile << std::fixed << std::setprecision(6) << v_bg_[k][y] << ' ';
		ofile << std::endl;
	}
}
