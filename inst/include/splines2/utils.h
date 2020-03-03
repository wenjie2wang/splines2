#ifndef SPLINE2_UTILS_H
#define SPLINE2_UTILS_H

#include <stdexcept>

#include "common.h"

namespace splines2 {

    // function to remove the first column of a matrix
    // would use const conference if possible:
    // see <https://github.com/RcppCore/Rcpp/issues/1045> for discussion
    inline rmat mat_wo_col1(rmat& x)
    {
        int x_ncol = x.ncol();
        if (x_ncol > 1) {
            rmat out { x(Rcpp::_, Rcpp::Range(1, x_ncol - 1)) };
            return out;
        }
        // else
        return rmat();
    }


}

#endif
