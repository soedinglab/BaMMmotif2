//
// Created by wanwan on 16.08.17.
//

#ifndef EM_H
#define EM_H

#include "BackgroundModel.h"
#include "MotifSet.h"

class EM {
    /**
     * This class is to optimize BaMM model using
     * expectation-maximization (EM) algorithm for
     * finding a local optimum
     */

public:

    EM( Motif* motif, BackgroundModel* bg, std::vector<Sequence*> seqs, float q );
    ~EM();

    int                     optimize();         // run EM optimization
    void                    print();            // print out optimized model v
    void					write( char* odir, std::string basename, bool ss );
                                                // write out the EM parameters such as n, pos, r

    void 					EStep();			// E-step
    void 					MStep();			// M-step

    float**                 getR();             // get the responsibility parameter r

private:

    Motif* 					motif_;				// motif to optimize within the EM
    BackgroundModel*		bg_;				// background model

    size_t 					K_;					// the order of the motif model
    size_t					W_;					// the width of the motif pattern
    float** 				A_;	        		// pseudo-count hyper-parameter for order k and motif position j
    size_t 					K_bg_;				// the order of the background model

    float** 				r_;		        	// responsibilities at all the positions in sequence n
                                                // Note: here the r_[n][0] indicates the responsibility of not having
                                                //      a motif on the sequence;
                                                //      r_[n][i] (for i > 0) indicates the responsibility of having a motif
                                                //      on position L-W+2-i
    float**					s_;					// log odds scores
    float*** 				n_;	            	// fractional counts n for (k+1)-mers y at motif position j
    float**					pos_;				// positional prior, pos[i][0] indicates the prior for no motif present on sequence i

    float 					q_; 				// hyper-parameter q specifies the fraction of sequences containing motif
    std::vector<Sequence*>	seqs_;				// copy positive sequences

    float 					llikelihood_ = 0.0f;// log likelihood for each iteration
    float					epsilon_ = 0.001f;	// threshold for likelihood convergence
    size_t					maxEMIterations_ = std::numeric_limits<size_t>::max();

    std::vector<size_t>		Y_;

};
#endif //EM_H
