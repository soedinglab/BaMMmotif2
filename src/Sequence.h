/*
 * Sequence.h
 *
 *  Created on: Apr 1, 2016
 *      Author: wanwan
 */

#ifndef SEQUENCE_H_
#define SEQUENCE_H_

class Sequence {

public:

	Sequence( char* sequence, int L, char* header ){
	    sequence_ = sequence;
	    L_ = L;
	    header_ = header;
	}
	~Sequence();

	unsigned int    getL(){ return L_; };				                        // get sequence length
	char*           getSequence(){ return sequence_; };	                        // get base sequence as alphabet encoding
	char* 	        getHeader(){ return header_; };		                        // get sequence header

	float 	        getIntensity(){ return intensity_; };                       // get measured sequence intensity
	float 	        getWeight(){ return weight_; };		                        // get sequence weight
	void 	        setIntensity( float intensity ){ intensity_ = intensity; };	// set measured sequence intensity
	void 	        setWeight( float weight ){ weight_ = weight; };				// set sequence weight

private:

	unsigned int    L_;			                        // sequence length
	char*		    sequence_;		                    // the base sequence as alphabet encoding
	char* 	        header_;			                // sequence header
	float 	        intensity_;		                    // measured sequence intensity from HT-SELEX data
	float 	        weight_;			                // sequence weight from HT-SELEX data, based on the intensity
};

#endif /* SEQUENCE_H_ */
