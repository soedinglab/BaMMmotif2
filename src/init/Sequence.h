#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <cstring>	// e.g. std::memcpy
#include <vector>

#include <stdint.h>	// e.g. uint8_t
#include <math.h>

#include "Alphabet.h"

class Sequence{

public:

	Sequence( uint8_t* sequence,
				size_t L,
				std::string header,
				bool singleStrand = false );
	~Sequence();

	uint8_t*		getSequence();
	size_t			getL();
	std::string		getHeader();


	float 	        getIntensity();
	float 	        getWeight();
	size_t*			getKmer();		    // get the value for 10-mer in the sequence

	void 	        setIntensity( float intensity );
	void 	        setWeight( float weight );

	void			print();		    // print out sequences

private:
					// append the sequence's reverse complement to the sequence
	void 			appendRevComp( uint8_t* sequence, size_t L );

	uint8_t*		sequence_;		    // sequence in alphabet encoding
	size_t			L_;				    // sequence length
	std::string		header_;		    // sequence header

	float			intensity_ = 0.0f;	// sequence intensity
	float			weight_ = 0.0f;	    // sequence weight calculated from its intensity
	size_t*			kmer_;
};

inline size_t* Sequence::getKmer(){
	return kmer_;
}

#endif /* SEQUENCE_H_ */
