#ifndef OBJFUN_H_
#define OBJFUN_H_

class ObjFun{

private:
    size_t          LW1_;
    size_t          posN_;
    float           N0_;
    float           beta2_;
    Eigen::VectorXf Ni_;
    Eigen::MatrixXf A_matrix_;

public:
    ObjFun( size_t LW1,
            size_t posN,
            float N0,
            float beta2,
            Eigen::VectorXf Ni,
            Eigen::MatrixXf A_matrix )
            :
            LW1_(LW1),
            posN_(posN),
            N0_(N0),
            beta2_(beta2),
            Ni_(Ni),
            A_matrix_(A_matrix)
    {};

    float operator()( Eigen::VectorXf& si, Eigen::VectorXf& grad ) {

        float Q;
        float p1 = 0.f;
        float p2 = 0.f;
        float sumexp = 0.f;

        // calculate objective function
        for (size_t i = 0; i < LW1_; i++) {
            sumexp += expf(si[i]);
        }
        for (size_t i = 0; i < LW1_; i++) {
            p1 += Ni_[i] * (si[i] - logf(sumexp));
        }
        for (size_t i = 2; i < LW1_; i++) {
            p2 += beta2_ * powf(si[i] - si[i - 1], 2.f);
        }

        Q = p2 - p1;
        //Q = p1 -p2;

        // update the gradients
        Eigen::VectorXf dot_product = A_matrix_ * si;
        for (size_t i = 0; i < LW1_; i++) {
            grad[i] = -(Ni_[i] - (posN_ - N0_) * expf(si[i]) / sumexp -
                        beta2_ * dot_product[i]);
        }

        return Q;
    }
};

#endif // end OBJFUN_H_