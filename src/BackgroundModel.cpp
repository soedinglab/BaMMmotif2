/*
 * BackgroundModel.cpp
 *
 *  Created on: Apr 18, 2016
 *      Author: administrator
 */

#include "BackgroundModel.h"
#include "SequenceSet.h"

BackgroundModel::BackgroundModel(){
	// allocate memory for v
	v_ = ( float** )calloc( ( Global::modelOrder+1 ), sizeof( float* ) );
	for( int k = 0; k <= Global::modelOrder; k++ ){
		v_[k] = ( float* )calloc( pow( Alphabet::getSize(), k+1 ), sizeof( float ) );
	}
}

BackgroundModel::~BackgroundModel(){
	if( v_ ) free( v_ );
}

void BackgroundModel::init(){				// the full negSeqSet is used

    // initialize n
    int** n;
	n = new int*[Global::modelOrder + 1];
	for( int k = 0; k < Global::modelOrder + 1; k++ ){
		int sizeY = int( pow( Alphabet::getSize(), k + 1 ) );
		n[k] = new int[sizeY];
		for( int y = 0; y < pow( Alphabet::getSize(), k + 1 ); y++ ){
			n[k][y] = 0;
		}
	}

	// count n
	for( unsigned int i = 0; i < Global::negSequenceSet->getN(); i++ ){
		Sequence sequence = Global::negSequenceSet->getSequences()[i];

		int l;								// cut seq by half when revcomp is true
		if( Global::revcomp )
			l = sequence.getL()/2;
		else
			l = sequence.getL();

		for( int k = 0; k < Global::modelOrder + 1; k++ ){

			for( int j = k; j < l; j++ ){
				int y = 0;
				for( int i = 0; i <= k; i++ )
					y += int( pow( Alphabet::getSize(), i ) * ( sequence.getSequence()[j-k+i] - 1 ) );
				n[k][y]++;
			}
		}
	}

	if( Global::verbose ){
		printf( " ___________________________\n"
				"|                           |\n"
				"| k-mer counts for bgModel  |\n"
				"|___________________________|\n\n" );
		for( int k = 0; k < Global::modelOrder + 1; k++ ){
			std::cout << "k = " << k << std::endl;
			for( int y = 0; y < pow( Alphabet::getSize(), k + 1 ); y++ ){
				std::cout << std::fixed << std::setprecision(4) << n[k][y] << '\t';
			}
			std::cout << std::endl;
		}
	}

    // calculate v from k-mer counts n
	calculateV( n );

	if( Global::verbose ) print();

	// free memory
	for( int i = 0; i < Global::modelOrder + 1; i++ ){
		delete[] n[i];
	}
	delete[] n;

}

void BackgroundModel::init( std::vector<int> folds ){		// for training sets in cross-validation

    // initialize n
    int** n;
	n = new int*[Global::modelOrder + 1];
	for( int k = 0; k < Global::modelOrder + 1; k++ ){
		int sizeY = int( pow( Alphabet::getSize(), k + 1 ) );
		n[k] = new int[sizeY];
		for( int y = 0; y < pow( Alphabet::getSize(), k + 1 ); y++ ){
			n[k][y] = 0;
		}
	}

	for( unsigned int f = 0; f < folds.size(); f++ ){
		// count n
		for( unsigned int r = 0; r < Global::negFoldIndices[f].size(); r++ ){
			unsigned int idx = Global::negFoldIndices[f][r];

			std::cout << "For negSequence " << idx << ":" << std::endl;

			Sequence sequence = Global::negSequenceSet->getSequences()[idx];

			int l;							// cut seq by half when revcomp is true
			if( Global::revcomp )
				l = sequence.getL()/2;
			else
				l = sequence.getL();

			for( int k = 0; k < Global::modelOrder + 1; k++ ){
				std::cout << "k = " << k << std::endl;
				for( int j = k; j < l; j++ ){
					int y = 0;
					for( int i = 0; i <= k; i++ )
						y += int( pow( Alphabet::getSize(), i ) * ( sequence.getSequence()[j-k+i] - 1 ) );
					n[k][y]++;
					std::cout << n[k][y] << '\t';
				}
				std::cout << std::endl;
			}
			std::cout << std::endl << std::endl;
		}
		std::cout << std::endl << std::endl;
	}

	if( Global::verbose ){
		printf( " ___________________________\n"
				"|                           |\n"
				"| k-mer counts for bgModel  |\n"
				"|___________________________|\n\n" );
		for( int k = 0; k < Global::modelOrder + 1; k++ ){
			std::cout << "k = " << k << std::endl;
			for( int y = 0; y < pow( Alphabet::getSize(), k + 1 ); y++ ){
				std::cout << std::fixed << std::setprecision(4) << n[k][y] << '\t';
			}
			std::cout << std::endl;
		}
	}

    // calculate v from k-mer counts n
	calculateV( n );

	if( Global::verbose ) print();

	// free memory
	for( int i = 0; i < Global::modelOrder + 1; i++ ){
		delete[] n[i];
	}
	delete[] n;
}

void BackgroundModel::calculateV( int** n ){

	int* sum = new int[Global::modelOrder + 1];
	for( int k = 0; k < Global::modelOrder + 1; k++ ){
		sum[k] = 0;
		for( int y = 0; y < pow( Alphabet::getSize(), k + 1 ); y++ ){
			sum[k] += n[k][y] ;
		}
	}

	/*
	for( int k = 0; k < Global::modelOrder + 1; k++ ){
		for( int y = 0; y < pow( Alphabet::getSize(), k + 1 ); y++ ){
			f_[k][y] = n[k][y] * ( k + 1 ) / float( sum[k] );
		}
	}
	*/

	// for k = 0:
	for( int y = 0, k = 0; y < pow( Alphabet::getSize(), k + 1 ); y++ ){
		v_[k][y] = n[k][y] * ( k + 1 ) / float( sum[k] );
	}

	// for k >0:
	for( int k = 1; k < Global::modelOrder + 1; k++ ){
		for( int y = 0; y < pow( Alphabet::getSize(), k + 1 ); y++ ){
			int y2 = y / pow( Alphabet::getSize(), k );
			int yk = y % int( pow( Alphabet::getSize(), k ) );
			v_[k][y] = ( n[k][y] + Global::modelAlpha[k] * v_[k-1][y2] ) / ( n[k-1][yk] + Global::modelAlpha[k] );
		}
	}
	delete[] sum;
}

float** BackgroundModel::getV(){
    return v_;
}

void BackgroundModel::print(){

	printf( " ___________________________\n"
			"|                           |\n"
			"| probabilities for bgModel |\n"
			"|___________________________|\n\n" );
	for( int k = 0; k < Global::modelOrder + 1; k++ ){
		std::cout << "when k = " << k << std::endl;
		for( int y = 0; y < pow( Alphabet::getSize(), k + 1 ); y++ ){
			std::cout << std::fixed << std::setprecision(4) << v_[k][y] << '\t';
		}
		std::cout << std::endl;
	}
}

void BackgroundModel::write(){
	std::string opath = std::string( Global::outputDirectory )  + '/' + std::string( Global::posSequenceBasename ) + ".probsBg";
	std::ofstream ofile( opath.c_str() );
	for( int k = 0; k < Global::modelOrder + 1; k++ ){
		for( int y = 0; y < pow( Alphabet::getSize(), k+1 ); y++ ){
			ofile << std::fixed << std::setprecision(4) << v_[k][y] << '\t';
		}
		ofile << std::endl;
	}
}
