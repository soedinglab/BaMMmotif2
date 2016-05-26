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
	// allocate memory for v
	v_ = ( float* )calloc( Y_, sizeof( float ) );
	std::cout << "Number of (k+1)-mers: " << Y_ << std::endl;
}

BackgroundModel::~BackgroundModel(){
	if( v_ ) free( v_ );
}

void BackgroundModel::init( std::vector<int> folds = std::vector<int> () ){
    // initialize n_occ[y]
    int* n_occ;
	n_occ = new int[Y_];
	for( int y = 0; y < Y_; y++ )
		n_occ[y] = 0;

	// count n_occ[y]
	if( folds.size() == 0 ){
		// when full negSeqSet is applied
		for( unsigned int n = 0; n < Global::negSequenceSet->getN(); n++ ){
			Sequence seq = Global::negSequenceSet->getSequences()[n];
			for( int k = 0; k < Global::bgModelOrder + 1; k++ ){
				for( int i = k; i < seq.getL(); i++ ){
					int y = seq.extractKmerbg( i, k );
					if( y >= 0 )				// skip the non-set alphabets, such as N
						n_occ[y]++;
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
							n_occ[y]++;
					}
				}
			}
		}
	}

	if( Global::verbose ){		// Only for testing
		printf( " ___________________________\n"
				"|                           |\n"
				"| k-mer counts for bgModel  |\n"
				"|___________________________|\n\n" );

		for( int y = 0; y < Y_; y++ )
			std::cout << n_occ[y] << '\t';
		std::cout << std::endl;
	}

    // calculate v from k-mer counts n_occ
	calculateV( n_occ );

	if( Global::verbose ) print();

	// free memory
	delete[] n_occ;
}

void BackgroundModel::calculateV( int* n_occ ){

	// sum of the counts n[y]
	float sum = 0.0f;
	for( int y = 0; y < powA_[1]; y++ )
		sum += float( n_occ[y] );

	// for k = 0: v = freqs
	for( int y = 0; y < powA_[1]; y++ )
		v_[y] = float( n_occ[y] ) / sum ;

	// for k > 0:
	for( int y = powA_[1]; y < Y_; y++ ){
		int y2 = y / powA_[1] - 1;						// cut off the first nucleotide, e.g. from ACGT to CGT
		for( int k = 1; k < Global::bgModelOrder; k++ ){
			int yk = y % powA_[k];						// cut off the last nucleotide, e.g. from ACGT to ACG
			v_[y] = ( n_occ[y] + Global::bgModelAlpha * v_[y2] ) / ( n_occ[yk] + Global::bgModelAlpha );
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
		for( int i = 0; i < powA_[k+1]; i ++, y++ )
			std::cout << std::fixed << std::setprecision(4) << v_[y] << '\t';
		std::cout << std::endl;
	}
}

void BackgroundModel::write(){
	std::string opath = std::string( Global::outputDirectory )  + '/' + std::string( Global::posSequenceBasename ) + ".probsBg";
	std::ofstream ofile( opath.c_str() );
	for( int k = 0, y = 0; k < Global::bgModelOrder + 1; k++){
		for( int i = 0; i < powA_[k+1]; i ++, y++ )
			ofile << std::fixed << std::setprecision(4) << v_[y] << '\t';
		ofile << std::endl;
	}
}
