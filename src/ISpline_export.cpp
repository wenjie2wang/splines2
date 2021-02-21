//
// R package splines2 by Wenjie Wang and Jun Yan
// Copyright (C) 2016-2021
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
Rcpp::NumericMatrix rcpp_iSpline(
    const arma::vec& x,
    const unsigned int df,
    const unsigned int degree,
    const arma::vec& internal_knots,
    const arma::vec& boundary_knots,
    const unsigned int derivs,
    const bool integral,
    const bool complete_basis = true
    )
{
    splines2::ISpline is_obj;
    Rcpp::NumericMatrix out;
    // if df > 0 and knots are not specified
    // auto set internal knots based on df
    if (df > 0 && internal_knots.n_elem == 0) {
        // compute actual spline degree of freedom
        const unsigned int wo_intercept {
            static_cast<unsigned int>(! complete_basis)
        };
        unsigned int spline_df { df + wo_intercept };
        is_obj = splines2::ISpline(x, spline_df, degree, boundary_knots);
    } else {
        // else ignore df
        is_obj = splines2::ISpline(x, internal_knots, degree, boundary_knots);
    }
    if (integral && derivs == 0) {
        // integrals
        out = splines2::arma2rmat(is_obj.integral(complete_basis));
    } else if ((! integral && derivs == 0) || (integral && derivs == 1)) {
        // basis functions
        out = splines2::arma2rmat(is_obj.basis(complete_basis));
    } else {
        // derivatives
        out = splines2::arma2rmat(
            is_obj.derivative(derivs - static_cast<unsigned int>(integral),
                              complete_basis)
            );
    }
    // add attributes
    out.attr("dimnames") = Rcpp::List::create(
        R_NilValue, splines2::char_seq_len(out.ncol())
        );
    out.attr("x") = splines2::arma2rvec(x);
    out.attr("degree") = static_cast<int>(is_obj.get_degree());
    out.attr("knots") = splines2::arma2rvec(is_obj.get_internal_knots());
    out.attr("Boundary.knots") =
        splines2::arma2rvec(is_obj.get_boundary_knots());
    out.attr("derivs") = static_cast<int>(derivs);
    out.attr("intercept") = complete_basis;
    return out;
}
