//
// Created by wanwan on 18.09.17.
//

#ifndef GSIMU_H_
#define GSIMU_H_
#include <string.h>	                            // e.g. strdup

#include "../init/Alphabet.h"
#include "../init/SequenceSet.h"

class GSimu{

public:

    static char*		outputDirectory; 		// output directory

    // sequence set options
    static char*		sequenceFilename;	    // filename of input sequence FASTA file
    static std::string	sequenceBasename;	    // basename of input sequence FASTA file
    static SequenceSet*	sequenceSet;			// input sequence set
    static float		q;						// prior probability for a positive sequence to contain a motif
    static char* 		alphabetType;			// provide alphabet type

    // initial model(s) options
    static char*		initialModelFilename;	// filename of initial model
    static std::string	initialModelTag;		// tag for initializing the model
    static std::string  fileExtension;          // extended filename for output

    // model options
    static size_t		modelOrder;				// model order
    static std::vector<float> modelAlpha;		// initial alphas
    static float		modelBeta;				// alpha_k = beta x gamma^k for k > 0
    static float		modelGamma;
    static std::vector<size_t>	addColumns;		// add columns to the left and right of init used to initialize Markov model

    // background model options
    static size_t		bgModelOrder;			// background model order, defaults to 2
    static std::vector<float> bgModelAlpha;		// background model alpha

    // simulation options
    static size_t		mFold;					// number of sampled sequences as multiple of input sequences
    static size_t		sOrder;					// k-mer order for sampling negative sequence set
    static bool         sampleBgset;
    static bool         maskSeqset;
    static bool         embedSeqset;
    static size_t       at;                     // position for embedding the motif

    static void         init( int nargs, char* args[] );
    static void         destruct();

private:

    static int	        readArguments( int nargs, char* args[] );
    static void	        printHelp();
};

#endif //GSIMU_H_
