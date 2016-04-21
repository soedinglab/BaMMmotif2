/*
 * Sequence.h
 *
 *  Created on: Apr 1, 2016
 *      Author: wanwan
 */

#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <stdint.h>		// uint8_t

class Sequence {

public:

	Sequence( uint8_t* sequence, int L, char* header );
	~Sequence();

	unsigned int    getL();				                // get sequence length
	uint8_t*        getSequence();	                   	// get base sequence as alphabet encoding
	char* 	        getHeader();		                // get sequence header

	float 	        getIntensity();                     // get measured sequence intensity
	float 	        getWeight();		                // get sequence weight
	void 	        setIntensity( float intensity );	// set measured sequence intensity
	void 	        setWeight( float weight );			// set sequence weight

private:

	unsigned int    L_;			                        // sequence length
	uint8_t*	    sequence_;		                    // the base sequence as alphabet encoding
	char* 	        header_;			                // sequence header
	float 	        intensity_;		                    // measured sequence intensity from HT-SELEX data
	float 	        weight_;			                // sequence weight from HT-SELEX data, based on the intensity

	uint8_t* 		createRevComp( uint8_t* seq );		// create reverse complement sequences for each sequence
};

#endif /* SEQUENCE_H_ */
