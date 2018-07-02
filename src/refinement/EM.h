//
// Created by wanwan on 16.08.17.
//

#ifndef EM_H_
#define EM_H_

#include "../init/BackgroundModel.h"
#include "../init/MotifSet.h"

class EM {
    /**
     * This class is to optimize BaMM model using
     * expectation-maximization (EM) algorithm for
     * finding a local optimum
     */

public:

    EM( Motif* motif, BackgroundModel* bgModel, std::vector<Sequence*> seqs,
        bool optimizeQ = true, bool verbose = false, float f = 0.2f );
    ~EM();

    int                     optimize();         // run EM optimization
    int                     mask();             // improve the EM optimization
    void                    print();            // print out optimized model v
    void					write( char* odir, std::string basename, bool ss );
                                                // write out the EM parameters such as n, pos, r

    void 					EStep();			// E-step
    void 					MStep();			// M-step
    void                    optimizeQ();        // optimize the hyper-parameter q
    void                    initializePos();    // initialize positional prior pos_i


    float**                 getR();             // get the responsibility parameter r
    float                   getQ();             // get the optimized positional prior q
    void                    printR();           // print out the responsibility parameter r

private:

    Motif* 					motif_;				// motif to optimize within the EM
    BackgroundModel*		bgModel_;			// background model

    size_t 					K_;					// the order of the motif model
    size_t					W_;					// the width of the motif pattern
    float** 				A_;	        		// pseudo-count hyper-parameter for order k and motif position j
    size_t 					K_bg_;				// the order of the background model

    float** 				r_;		        	// responsibilities at all the positions in sequence n
                                                // Note: here the r_[n][0] indicates the responsibility of not having
                                                //      a motif on the n'th sequence;
                                                //      r_[n][i] (for i > 0) indicates the responsibility of having a motif
                                                //      at position L-W+2-i on the n'th sequence
    float**					s_;					// log odds scores
    float*** 				n_;	            	// fractional counts n for (k+1)-mers y at motif position j
    float**					pos_;				// positional prior, pos[i][0] indicates the prior for no motif present on sequence i

    float 					q_; 				// hyper-parameter q specifies the fraction of sequences containing motif
    float                   f_;                 // fraction of sequences to be masked
    std::vector<Sequence*>	seqs_;				// copy positive sequences

    float 					llikelihood_        = 0.0f;     // log likelihood for each iteration
    float					epsilon_            = 0.01f;	// threshold for parameter v convergence
    size_t					maxEMIterations_    = 1000;
    bool                    optimizeQ_;

    bool                    verbose_;           // show the output of each EM iteration
    std::vector<size_t>		Y_;

};

#endif //EM_H_
