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
Rcpp::NumericMatrix rcpp_mSpline_basis(
    const arma::vec& x,
    const unsigned int df,
    const unsigned int degree,
    const arma::vec& internal_knots,
    const arma::vec& boundary_knots,
    const bool complete_basis = true
    )
{
    const unsigned int wo_intercept {
        static_cast<unsigned int>(! complete_basis)
    };
    splines2::MSpline ms_obj;
    // if df > 0 and knots are not specified
    // auto set internal knots based on df
    if (df > 0 && internal_knots.n_elem == 0) {
        // compute actual spline degree of freedom
        unsigned int spline_df { df + wo_intercept };
        ms_obj = splines2::MSpline(x, spline_df, degree, boundary_knots);
    } else {
        // else ignore df
        ms_obj = splines2::MSpline(x, internal_knots, degree, boundary_knots);
    }
    Rcpp::NumericMatrix out {
        splines2::arma2rmat(ms_obj.basis(complete_basis))
            };
    // add attributes
    out.attr("dimnames") = Rcpp::List::create(
        R_NilValue, splines2::char_seq_len(out.ncol())
        );
    out.attr("x") = splines2::arma2rvec(x);
    out.attr("degree") = static_cast<int>(ms_obj.get_degree());
    out.attr("knots") = splines2::arma2rvec(ms_obj.get_internal_knots());
    out.attr("Boundary.knots") =
        splines2::arma2rvec(ms_obj.get_boundary_knots());
    out.attr("intercept") = complete_basis;
    return out;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_mSpline_derivative(
    const arma::vec& x,
    const unsigned int derivs,
    const unsigned int df,
    const unsigned int degree,
    const arma::vec& internal_knots,
    const arma::vec& boundary_knots,
    const bool complete_basis = true
    )
{
    const unsigned int wo_intercept {
        static_cast<unsigned int>(! complete_basis)
    };
    splines2::MSpline ms_obj;
    // if df > 0 and knots are not specified
    // auto set internal knots based on df
    if (df > 0 && internal_knots.n_elem == 0) {
        // compute actual spline degree of freedom
        unsigned int spline_df { df + wo_intercept };
        ms_obj = splines2::MSpline(x, spline_df, degree, boundary_knots);
    } else {
        // else ignore df
        ms_obj = splines2::MSpline(x, internal_knots, degree, boundary_knots);
    }
    Rcpp::NumericMatrix out {
        splines2::arma2rmat(ms_obj.derivative(derivs, complete_basis))
            };
    // add attributes
    out.attr("dimnames") = Rcpp::List::create(
        R_NilValue, splines2::char_seq_len(out.ncol())
        );
    out.attr("x") = splines2::arma2rvec(x);
    out.attr("derivs") = static_cast<int>(derivs);
    out.attr("degree") = static_cast<int>(ms_obj.get_degree());
    out.attr("knots") = splines2::arma2rvec(ms_obj.get_internal_knots());
    out.attr("Boundary.knots") =
        splines2::arma2rvec(ms_obj.get_boundary_knots());
    out.attr("intercept") = complete_basis;
    return out;
}
