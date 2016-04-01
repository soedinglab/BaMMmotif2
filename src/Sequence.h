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

	Sequence( int* sequence, int L, char* header );
	~Sequence();

	int		getL();			// get sequence length
	int* 	getSequence();	// get base sequence as alphabet encoding
	char* 	getHeader();	// get sequence header
	float 	getIntensity;	// get measured sequence intensity
	float 	getWeight;		// get sequence weight

	float 	setIntensity;	// set measured sequence intensity
	float 	setWeight;		// set sequence weight

private:

	int 	L;				// sequence length
	char*	sequence;		// the base sequence as alphabet encoding
	char* 	header;			// sequence header
	float 	intensity;		// measured sequence intensity
	float 	weight;			// sequence weight
};



#endif /* SEQUENCE_H_ */
