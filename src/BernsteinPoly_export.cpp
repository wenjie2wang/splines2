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

#include <RcppArmadillo.h>
// [[Rcpp::plugins(cpp11)]]

#include <splines2Armadillo.h>

// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_bernsteinPoly(
    const arma::vec& x,
    const unsigned int degree,
    const unsigned int derivs = 0,
    const bool integral = false,
    const bool complete_basis = true
    )
{
    splines2::BernsteinPoly bp_obj { x, degree };
    Rcpp::NumericMatrix out;
    if (integral && derivs == 0) {
        // integrals
        out = splines2::arma2rmat(bp_obj.integral(complete_basis));
    } else if ((! integral && derivs == 0) || (integral && derivs == 1)) {
        // bases
        out = splines2::arma2rmat(bp_obj.basis(complete_basis));
    } else {
        // derivatives
        out = splines2::arma2rmat(
            bp_obj.derivative(derivs - static_cast<unsigned int>(integral),
                              complete_basis)
            );
    }
    // add attributes
    out.attr("dimnames") = Rcpp::List::create(
        R_NilValue, splines2::char_seq_len(out.ncol())
        );
    out.attr("x") = splines2::arma2rvec(x);
    out.attr("degree") = bp_obj.get_degree();
    out.attr("derivs") = derivs;
    out.attr("integral") = integral;
    out.attr("intercept") = complete_basis;
    return out;
}
