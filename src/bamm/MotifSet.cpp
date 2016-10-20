#include <fstream>		// std::fstream

#include "MotifSet.h"

MotifSet::MotifSet(){

	if( Global::BaMMpatternFilename != NULL ){

	    // scan file and conduct for IUPAC pattern p
		// * Motif* motif = new Motif( length(p) )
		// * motif.initFromIUPACPattern( p )
		// motifs.push_back( motif )

	} else if( Global::bindingSiteFilename != NULL ){

		std::ifstream file;
		file.open( Global::bindingSiteFilename, std::ifstream::in );
		std::string seq;
		int length;

		if( !file.good() ){
			std::cout << "Error: Cannot open bindingSitesFile sequence file: "
					<< Global::bindingSiteFilename << std::endl;
			exit( -1 );
		} else {
			getline( file, seq );					// get length of the first sequence
			length = ( int ) seq.length();
		}

		length += Global::addColumns.at( 0 ) + Global::addColumns.at( 1 );

		Motif* motif = new Motif( length );

		motif->initFromBindingSites( Global::bindingSiteFilename );

		motifs_.push_back( motif );

		N_ = 1;

	} else if( Global::PWMFilename != NULL ){

		// read file to calculate motif length
		// Motif* motif = new Motif( length )
		// motif.initFromPWM( file )
		// motifs.push_back( motif )

	} else if( Global::BaMMFilename != NULL ){

		// read file to calculate motif length
		// Motif* motif = new Motif( length )
		// motif.initFromBMM( file )
		// motifs.push_back( motif )
	}

//	if( Global::verbose ) print();
}

MotifSet::~MotifSet(){
    for( size_t i = 0; i < motifs_.size(); i++ ){
        delete motifs_[i];
    }
}

std::vector<Motif*> MotifSet::getMotifs(){
    return motifs_;
}

int MotifSet::getN(){
    return N_;
}

void MotifSet::print(){
	printf( " ____________________________\n"
			"|                            |\n"
			"| PROBABILITIES for MotifSet |\n"
			"|____________________________|\n\n" );

	for( int i = 0; i < N_; i++ ){
		printf( " ________________________________________\n"
				"|                                        |\n"
				"| INITIALIZED PROBABILITIES for Motif %d  |\n"
				"|________________________________________|\n\n", i+1 );
		motifs_[i]->print();
	}
}
void MotifSet::write(){

	/**
	 * save all the motifs learned by BaMM in two flat files:
	 * (1) posSequenceBasename.ihbcp: 		conditional probabilities after EM
	 * (2) posSequenceBasename.ihbp: 		probabilities of PWM after EM
	 */

	std::string opath = std::string( Global::outputDirectory )  + '/'
			+ std::string( Global::posSequenceBasename );

	// output conditional probabilities v[k][y][j] and probabilities prob[k][y][j]
	std::string opath_v = opath + ".ihbcp"; 	// inhomogeneous bamm conditional probabilities
	std::string opath_p = opath + ".ihbp";		// inhomogeneous bamm probabilities
	std::ofstream ofile_v( opath_v.c_str() );
	std::ofstream ofile_p( opath_p.c_str() );
	int i, j, k, y;

	std::vector<int> Y;
	for( k = 0; k < Global::modelOrder+2; k++ ){
		Y.push_back( ipow( Alphabet::getSize(), k ) );
	}

	for( i = 0; i < N_; i++ ){

//		ofile_v << "> motif " << i+1 <<":" << std::endl;
//		ofile_p << "> motif " << i+1 <<":" << std::endl;

		for( j = 0; j < motifs_[i]->getW(); j++ ){
			for( k = 0; k < Global::modelOrder+1; k++ ){
				for( y = 0; y < Y[k+1]; y++ ){
					ofile_v << std::scientific << std::setprecision(8) << motifs_[i]->getV()[k][y][j] << ' ';
					ofile_p << std::scientific << std::setprecision(8) << motifs_[i]->getP()[k][y][j] << ' ';
				}
				ofile_v << std::endl;
				ofile_p << std::endl;
			}
			ofile_v << std::endl;
			ofile_p << std::endl;
		}
	}
}
