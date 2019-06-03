//
// Created by wanwan on 16.08.17.
//

#ifndef EM_H_
#define EM_H_

#include <random>
#include <chrono>

#include "../init/BackgroundModel.h"
#include "../init/MotifSet.h"
#include "../Global/Global.h"

// todo: for positional prior
#include "../eigen3/Eigen/Dense"   // e.g. conjugate gradient solver
#include "../eigen3/Eigen/IterativeLinearSolvers"
#include "../eigen3/Eigen/Core"    // e.g. used for LBFGS++
#include "ObjFun.h"
#include "../LBFGS/LBFGS.h"

class EM {
    /**
     * This class is to optimize BaMM model using
     * expectation-maximization (EM) algorithm for
     * finding a local optimum
     */

public:

    EM( Motif* motif, BackgroundModel* bgModel,
        std::vector<Sequence*> seqs );
    ~EM();

    int                     optimize();         // run EM optimization
    int                     mask();             // improve the EM optimization
    void                    print();            // print out optimized model v
    void					write( char* odir,
                                   std::string basename );
                                                // write out the EM parameters such as n, pos, r

    void 					EStep();			// E-step
    void 					MStep();			// M-step
    void 					EStep_slow();		// E-step: slow version
    void 					MStep_slow();		// M-step: slow version

    void                    optimizeQ();        // optimize the hyper-parameter q

    void                    updatePrior();        // initialize the positional prior pos_i

    // todo: for positional prior
    void                    optimizePos();      // optimize the positional prior pos_i
    Eigen::MatrixXf         getAmatrix( size_t w );
    Eigen::MatrixXf         getBmatrix( size_t w );
    float                   obj_fun( Eigen::VectorXf& si, Eigen::VectorXf& grad );
    float*                  getPi();            // get positional prior pi

    float**                 getR();             // get the responsibility parameter r
    float                   getQ();             // get the optimized positional prior q
    void                    printR();           // print out the responsibility parameter r

private:

    Motif* 					motif_;				// motif to optimize within the EM
    BackgroundModel*		bgModel_;			// background model
    std::vector<Sequence*>	seqs_;				// copy positive sequences

    size_t 					K_;					// the order of the motif model
    size_t					W_;					// the width of the motif pattern
    float** 				A_;	        		// pseudo-count hyper-parameter for order k and motif position j

    float** 				r_;		        	// responsibilities at all the positions in sequence n
                                                // Note: here the r_[n][0] indicates the responsibility of not having
                                                //      a motif on the n'th sequence;
                                                //      r_[n][i] (for i > 0) indicates the responsibility of having a motif
                                                //      at position L-W+2-i on the n'th sequence
    float**					s_;					// log odds scores
    float*** 				n_;	            	// fractional counts n for (k+1)-mers y at motif position j
    float**					prior_;				// positional prior, pos[i][0] indicates the prior for no motif present on sequence i

    float 					q_; 				// hyper-parameter q specifies the fraction of sequences containing motif
    size_t                  padding_;
    float 					llikelihood_= 0.0f; // log likelihood for each iteration

    // todo: for positional prior
    float                   beta1_;             // hyper-parameter for smoothness on positional prior
    float                   beta2_;             // hyper-parameter for smoothness on positional prior
    float                   norm_;
    float*                  pi_;                // probability of a motif to start at position i on the longest sequence
    float                   N0_;

    Eigen::VectorXf         b_vector_;
    Eigen::VectorXf         si_;
    Eigen::VectorXf         Ni_;
    Eigen::MatrixXf         A_matrix_;

};

#endif //EM_H_
