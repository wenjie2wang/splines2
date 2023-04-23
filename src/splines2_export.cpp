//
// R package splines2 by Wenjie Wang and Jun Yan
// Copyright (C) 2016-2023
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
#include <splines2Armadillo.h>

template <typename T>
Rcpp::NumericMatrix bm_spline(
    const arma::vec& x,
    const unsigned int df,
    const unsigned int degree,
    const arma::vec& internal_knots,
    const arma::vec& boundary_knots,
    const bool complete_basis = true,
    const unsigned int derivs = 0,
    const bool integral = false
    )
{
    Rcpp::NumericMatrix out;
    T sp_obj;
    // if df > 0 and knots are not specified
    // auto set internal knots based on df
    if (df > 0 && internal_knots.n_elem == 0) {
        // compute actual spline degree of freedom
        const unsigned int wo_intercept {
            static_cast<unsigned int>(! complete_basis)
        };
        unsigned int spline_df { df + wo_intercept };
        sp_obj = T(x, spline_df, degree, boundary_knots);
    } else {
        // else ignore df
        sp_obj = T(x, internal_knots,
                                   degree, boundary_knots);
    }
    // 1) basis, 2) derivative, or 3) integral
    if (integral && derivs == 0) {
        // integrals
        out = splines2::arma2rmat(sp_obj.integral(complete_basis));
    } else if ((! integral && derivs == 0) || (integral && derivs == 1)) {
        // basis functions
        out = splines2::arma2rmat(sp_obj.basis(complete_basis));
    } else {
        // derivatives
        out = splines2::arma2rmat(
            sp_obj.derivative(derivs - static_cast<unsigned int>(integral),
                              complete_basis)
            );
    }
    // add attributes
    out.attr("dimnames") = Rcpp::List::create(
        R_NilValue, splines2::char_seq_len(out.ncol())
        );
    out.attr("x") = splines2::arma2rvec(x);
    out.attr("degree") = static_cast<int>(sp_obj.get_degree());
    out.attr("knots") = splines2::arma2rvec(sp_obj.get_internal_knots());
    out.attr("Boundary.knots") =
        splines2::arma2rvec(sp_obj.get_boundary_knots());
    out.attr("intercept") = complete_basis;
    out.attr("derivs") = static_cast<int>(derivs);
    out.attr("integral") = integral;
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_bSpline(
    const arma::vec& x,
    const unsigned int df,
    const unsigned int degree,
    const arma::vec& internal_knots,
    const arma::vec& boundary_knots,
    const bool complete_basis = true,
    const bool periodic = false,
    const unsigned int derivs = 0,
    const bool integral = false
    )
{
    Rcpp::NumericMatrix out;
    if (periodic) {
        out = bm_spline<splines2::PeriodicBSpline>(
            x, df, degree, internal_knots, boundary_knots,
            complete_basis, derivs, integral);
    } else {
        out = bm_spline<splines2::BSpline>(
            x, df, degree, internal_knots, boundary_knots,
            complete_basis, derivs, integral);
    }
    out.attr("periodic") = periodic;
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_mSpline(
    const arma::vec& x,
    const unsigned int df,
    const unsigned int degree,
    const arma::vec& internal_knots,
    const arma::vec& boundary_knots,
    const bool complete_basis = true,
    const bool periodic = false,
    const unsigned int derivs = 0,
    const bool integral = false
    )
{
    Rcpp::NumericMatrix out;
    if (periodic) {
        out = bm_spline<splines2::PeriodicMSpline>(
            x, df, degree, internal_knots, boundary_knots,
            complete_basis, derivs, integral);
    } else {
        out = bm_spline<splines2::MSpline>(
            x, df, degree, internal_knots, boundary_knots,
            complete_basis, derivs, integral);
    }
    out.attr("periodic") = periodic;
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_bernsteinPoly(
    const arma::vec& x,
    const unsigned int degree,
    const arma::vec& boundary_knots,
    const bool complete_basis = true,
    const unsigned int derivs = 0,
    const bool integral = false
    )
{
    splines2::BernsteinPoly bp_obj { x, degree, boundary_knots };
    Rcpp::NumericMatrix out;
    if (integral && derivs == 0) {
        // integrals
        out = splines2::arma2rmat(bp_obj.integral(complete_basis));
    } else if ((! integral && derivs == 0) || (integral && derivs == 1)) {
        // basis functions
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
    out.attr("degree") = static_cast<int>(bp_obj.get_degree());
    out.attr("Boundary.knots") = splines2::arma2rvec(
        bp_obj.get_boundary_knots()
        );
    out.attr("intercept") = complete_basis;
    out.attr("derivs") = static_cast<int>(derivs);
    out.attr("integral") = integral;
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_iSpline(
    const arma::vec& x,
    const unsigned int df,
    const unsigned int degree,
    const arma::vec& internal_knots,
    const arma::vec& boundary_knots,
    const bool complete_basis = false,
    const unsigned int derivs = 0,
    const bool integral = false
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
    out.attr("intercept") = complete_basis;
    out.attr("derivs") = static_cast<int>(derivs);
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_cSpline(
    const arma::vec& x,
    const unsigned int df,
    const unsigned int degree,
    const arma::vec& internal_knots,
    const arma::vec& boundary_knots,
    const bool complete_basis = false,
    const unsigned int derivs = 0
    )
{
    splines2::CSpline cs_obj;
    Rcpp::NumericMatrix out;
    // if df > 0 and knots are not specified
    // auto set internal knots based on df
    if (df > 0 && internal_knots.n_elem == 0) {
        // compute actual spline degree of freedom
        const unsigned int wo_intercept {
            static_cast<unsigned int>(! complete_basis)
        };
        unsigned int spline_df { df + wo_intercept };
        cs_obj = splines2::CSpline(x, spline_df, degree, boundary_knots);
    } else {
        // else ignore df
        cs_obj = splines2::CSpline(x, internal_knots, degree, boundary_knots);
    }
    if (derivs > 0) {
        out = splines2::arma2rmat(cs_obj.derivative(derivs, complete_basis));
    } else {
        out = splines2::arma2rmat(cs_obj.basis(complete_basis));
    }
    // add attributes
    out.attr("dimnames") = Rcpp::List::create(
        R_NilValue, splines2::char_seq_len(out.ncol())
        );
    out.attr("x") = splines2::arma2rvec(x);
    out.attr("degree") = static_cast<int>(cs_obj.get_degree());
    out.attr("knots") = splines2::arma2rvec(cs_obj.get_internal_knots());
    out.attr("Boundary.knots") =
        splines2::arma2rvec(cs_obj.get_boundary_knots());
    out.attr("intercept") = complete_basis;
    out.attr("derivs") = static_cast<int>(derivs);
    out.attr("scales") = splines2::arma2rvec(cs_obj.get_scales());
    return out;
}


template <typename T>
Rcpp::NumericMatrix template_naturalSpline(
    const arma::vec& x,
    const unsigned int df,
    const arma::vec& internal_knots,
    const arma::vec& boundary_knots,
    const bool complete_basis = false,
    const unsigned int derivs = 0,
    const bool integral = false
    )
{
    T ns_obj;
    // if df > 0 and knots are not specified
    // auto set internal knots based on df
    if (df > 0 && internal_knots.n_elem == 0) {
        // compute actual spline degree of freedom
        const unsigned int wo_intercept {
            static_cast<unsigned int>(! complete_basis)
        };
        unsigned int spline_df { df + wo_intercept };
        ns_obj = T(x, spline_df, boundary_knots);
    } else {
        // else ignore df
        ns_obj = T(x, internal_knots, boundary_knots);
    }
    Rcpp::NumericMatrix out;
    // 1) basis, 2) derivative, or 3) integral
    if (integral && derivs == 0) {
        // integrals
        out = splines2::arma2rmat(ns_obj.integral(complete_basis));
    } else if ((! integral && derivs == 0) || (integral && derivs == 1)) {
        // basis functions
        out = splines2::arma2rmat(ns_obj.basis(complete_basis));
    } else {
        // derivatives
        out = splines2::arma2rmat(
            ns_obj.derivative(derivs - static_cast<unsigned int>(integral),
                              complete_basis)
            );
    }
    // add attributes
    out.attr("dimnames") = Rcpp::List::create(
        R_NilValue, splines2::char_seq_len(out.ncol())
        );
    out.attr("x") = splines2::arma2rvec(x);
    out.attr("knots") = splines2::arma2rvec(ns_obj.get_internal_knots());
    out.attr("Boundary.knots") =
        splines2::arma2rvec(ns_obj.get_boundary_knots());
    out.attr("intercept") = complete_basis;
    out.attr("derivs") = static_cast<int>(derivs);
    out.attr("integral") = integral;
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_naturalSpline(
    const arma::vec& x,
    const unsigned int df,
    const arma::vec& internal_knots,
    const arma::vec& boundary_knots,
    const bool complete_basis = false,
    const unsigned int derivs = 0,
    const bool integral = false
    )
{
    return template_naturalSpline<splines2::NaturalSpline>(
        x, df, internal_knots, boundary_knots,
        complete_basis, derivs, integral
        );
}

// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_nsk(
    const arma::vec& x,
    const unsigned int df,
    const arma::vec& internal_knots,
    const arma::vec& boundary_knots,
    const bool complete_basis = false,
    const unsigned int derivs = 0,
    const bool integral = false
    )
{
    return template_naturalSpline<splines2::NaturalSplineK>(
        x, df, internal_knots, boundary_knots,
        complete_basis, derivs, integral
        );
}
