/*
 * FDR.h
 *
 *  Created on: Apr 1, 2016
 *      Author: wanwan
 */


#ifndef FDR_H_
#define FDR_H_

#include "../shared/BackgroundModel.h"
#include "Motif.h"
#include "EM.h"

class FDR {

public:

	FDR( const Motif& motif );
	~FDR();

	void evaluateMotif();
	void print();
	void write();

private:

	const Motif&	motif_;				// motif learned on full SequenceSet

	std::vector< std::vector<int> >	testFoldIndices_;	// sequence indices for each cross-validation fold
	std::vector< std::vector<int> >	trainFoldIndices_;	// sequence indices for each cross-validation fold

	float*	 		posS_all_;			// log odds scores on positive test SequenceSets
	float* 			negS_all_;			// log odds scores on negative test SequenceSets
	float* 			posS_max_;			// largest log odds scores on positive test SequenceSets
	float* 			negS_max_;			// largest log odds scores on negative test SequenceSets

	float* 			P_ZOOPS_;			// precision for ZOOPS model
	float* 			R_ZOOPS_;			// recall for ZOOPS model

	float*			P_MOPS_;			// precision for MOPS model
	float*			R_MOPS_;			// recall for MOPS model
	float*			FP_MOPS_;			// false positive values for MOPS model
	float*			TFP_MOPS_;			// true and false positives value for MOPS model

	std::vector<Sequence>	posSeqs_ = Global::posSequenceSet->getSequences();
	int 			W_;					// motif length

					// generate background sample sequence set
	SequenceSet*	sampleSequenceSet( EM* em );
	Sequence*		sampleSequence();

					// generate fold indices for test and training sets
	void			generateFoldIndices( int posN, int fold );

					// score Global::posSequenceSet using Global::foldIndices and folds
	void 	        scoreSequenceSet( Motif* motif, BackgroundModel* bg, std::vector<int> folds );

					// score SequenceSet sequences
	void 		    scoreSequenceSet( Motif* motif, BackgroundModel* bg, SequenceSet* seqSet );

	void 		    calculatePR();
					// Quick sort algorithm in descending order
	void			quickSort( float* arr, int left, int right );
};

inline void quickSort( float* arr, int left, int right ) {

	int i = left, j = right;
	float tmp;
	float pivot = arr[( left + right ) / 2];

	/* partition */
	while( i <= j ){
		while( arr[i] - pivot > 0 )	i++;
		while( arr[j] - pivot < 0 )	j--;
		if( i <= j ){
			tmp = arr[i];
			arr[i] = arr[j];
			arr[j] = tmp;
			i++;
			j--;
		}
	}

	/* recursion */
	if( left < j )	quickSort( arr, left, j );
	if( i < right )	quickSort( arr, i, right );
}

#endif /* FDR_H_ */

