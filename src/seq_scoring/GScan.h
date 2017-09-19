//
// Created by wanwan on 04.09.17.
//

#ifndef GSCAN_H_
#define GSCAN_H_

#include <string.h>	                            // e.g. strdup

#include "../init/Alphabet.h"
#include "../init/SequenceSet.h"

class GScan{

public:

    static char*		outputDirectory; 		// output directory

    // sequence set options
    static char*		posSequenceFilename;	// filename of positive sequence FASTA file
    static std::string	posSequenceBasename;	// basename of positive sequence FASTA file
    static SequenceSet*	posSequenceSet;			// positive sequence set
    static float		q;						// prior probability for a positive sequence to contain a motif
    static bool			ss;						// only search on single strand sequences
    static char*		negSequenceFilename;	// filename of negative sequence FASTA file
    static SequenceSet*	negSequenceSet;			// negative sequence set
    static char* 		alphabetType;			// provide alphabet type

    // initial model(s) options
    static char*		initialModelFilename;	// filename of initial model
    static std::string	initialModelTag;		// tag for initializing the model
    static size_t		num;					// number of init that are to be optimized

    // model options
    static size_t		modelOrder;				// model order
    static std::vector<float> modelAlpha;		// initial alphas
    static float		modelBeta;				// alpha_k = beta x gamma^k for k > 0
    static float		modelGamma;
    static std::vector<size_t>	addColumns;		// add columns to the left and right of init used to initialize Markov init
    static bool			interpolateBG;			// calculate prior probabilities from lower-order probabilities
                                                // instead of background frequencies of mono-nucleotides

    // background model options
    static size_t		bgModelOrder;			// background model order, defaults to 2
    static std::vector<float> bgModelAlpha;		// background model alpha

    static float        cutoff;                 // cutoff for scanning motifs

    static void         init( int nargs, char* args[] );
    static void         destruct();

private:

    static int	        readArguments( int nargs, char* args[] );
    static void	        printHelp();
};
#endif // GSCAN_H_
