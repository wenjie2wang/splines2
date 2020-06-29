#ifndef SPLINES2_COMMON_H
#define SPLINES2_COMMON_H

#include <RcppArmadillo.h>

namespace splines2arma {

    typedef arma::vec rvec;
    typedef arma::mat rmat;
    typedef arma::ivec ivec;
    typedef arma::uvec uvec;

    // nan vector
    inline rvec nan_vec() {
        rvec out { arma::zeros(1) };
        out.fill(arma::datum::nan);
        return out;
    }

}  // splines2arma

#endif /* SPLINES2_COMMON_H */
