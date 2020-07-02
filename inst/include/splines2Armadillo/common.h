//
// R package splines2 by Wenjie Wang and Jun Yan
// Copyright (C) 2016-2020
//
// This file is part of the R package splines2.
//
// The R package splines2 is free software: You can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or any later
// version (at your option). See the GNU General Public License at
// <https://www.gnu.org/licenses/> for details.
//
// The R package splines2 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//

#ifndef SPLINES2_COMMON_H
#define SPLINES2_COMMON_H

#include <RcppArmadillo.h>

namespace splines2 {

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

}  // splines2

#endif /* SPLINES2_COMMON_H */
