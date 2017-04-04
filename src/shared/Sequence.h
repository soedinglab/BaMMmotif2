#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <cstring>	// e.g. std::memcpy
#include <vector>

#include <stdint.h>	// e.g. uint8_t
#include <math.h>
#include "Alphabet.h"

class Sequence{

public:

	Sequence( uint8_t* sequence, int L, std::string header, std::vector<int> Y, bool revcomp = false );
	~Sequence();

	uint8_t*        getSequence();
	int    			getL();
	std::string		getHeader();


	float 	        getIntensity();
	float 	        getWeight();
	int*			getKmer();						// get the value for 8-mer in the sequence

	void 	        setIntensity( float intensity );
	void 	        setWeight( float weight );

	void			print();						// print sequences

private:

	void 			appendRevComp( uint8_t* sequence, int L ); // append the sequence's reverse complement to the sequence

	uint8_t*		sequence_;						// sequence in alphabet encoding
	int				L_;								// sequence length
	std::string		header_;						// sequence header

	float			intensity_ = 0.0f;				// sequence intensity
	float			weight_ = 0.0f;					// sequence weight calculated from its intensity
	int*			kmer_;							// get the value for 8-mers
	int* 			Y_;								// contains 1 at position 0
													// and the number of oligomers y for increasing order k at positions k+1
													// e.g.
													// alphabet size_ = 4: Y_ = 4^0 4^1 4^2 ... 4^15 < std::numeric_limits<int>::max()
													// limits the length of oligomers to 15 (and the order to 14)
};

inline int* Sequence::getKmer(){
	return kmer_;
}

#endif /* SEQUENCE_H_ */
