/*
 * Sequence.h
 *
 *  Created on: Apr 1, 2016
 *      Author: wanwan
 */

#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <string>		// std::string

#include <stdint.h>		// uint8_t

class Sequence {

public:

	Sequence( uint8_t* sequence, int L, std::string header );
	~Sequence();

	int    			getL();				                // get sequence length
	uint8_t*        getSequence();	                   	// get base sequence as alphabet encoding
	std::string		getHeader();		                // get sequence header
	float 	        getIntensity();                     // get measured sequence intensity
	float 	        getWeight();		                // get sequence weight
	void 	        setIntensity( float intensity );	// set measured sequence intensity
	void 	        setWeight( float weight );			// set sequence weight
	int				extractKmer( int i, int k );		// extract k-mer from positions(i-k, ..., i) of the sequence
	int				extractKmerbg( int i, int k );		// extract k-mer from positions(i-k, ..., i) of the sequence for bgModel

private:

	int    			L_;			                        // sequence length
	uint8_t*	    sequence_;		                    // the base sequence as alphabet encoding
	std::string		header_;			                // sequence header
	float 	        intensity_;		                    // measured sequence intensity from HT-SELEX data
	float 	        weight_;			                // sequence weight from HT-SELEX data, based on the intensity

	uint8_t* 		createRevComp( uint8_t* sequence, int L );	// create reverse complement sequences for each sequence
};

#endif /* SEQUENCE_H_ */
