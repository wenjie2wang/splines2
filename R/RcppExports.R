# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

rcpp_bSpline <- function(x, df, degree, internal_knots, boundary_knots, complete_basis = TRUE, periodic = FALSE, derivs = 0L, integral = FALSE) {
    .Call('_splines2_rcpp_bSpline', PACKAGE = 'splines2', x, df, degree, internal_knots, boundary_knots, complete_basis, periodic, derivs, integral)
}

rcpp_mSpline <- function(x, df, degree, internal_knots, boundary_knots, complete_basis = TRUE, periodic = FALSE, derivs = 0L, integral = FALSE) {
    .Call('_splines2_rcpp_mSpline', PACKAGE = 'splines2', x, df, degree, internal_knots, boundary_knots, complete_basis, periodic, derivs, integral)
}

rcpp_bernsteinPoly <- function(x, degree, boundary_knots, complete_basis = TRUE, derivs = 0L, integral = FALSE) {
    .Call('_splines2_rcpp_bernsteinPoly', PACKAGE = 'splines2', x, degree, boundary_knots, complete_basis, derivs, integral)
}

rcpp_iSpline <- function(x, df, degree, internal_knots, boundary_knots, complete_basis = FALSE, derivs = 0L, integral = FALSE) {
    .Call('_splines2_rcpp_iSpline', PACKAGE = 'splines2', x, df, degree, internal_knots, boundary_knots, complete_basis, derivs, integral)
}

rcpp_cSpline <- function(x, df, degree, internal_knots, boundary_knots, complete_basis = FALSE, derivs = 0L) {
    .Call('_splines2_rcpp_cSpline', PACKAGE = 'splines2', x, df, degree, internal_knots, boundary_knots, complete_basis, derivs)
}

rcpp_naturalSpline <- function(x, df, internal_knots, boundary_knots, trim = 0.0, complete_basis = FALSE, derivs = 0L, integral = FALSE) {
    .Call('_splines2_rcpp_naturalSpline', PACKAGE = 'splines2', x, df, internal_knots, boundary_knots, trim, complete_basis, derivs, integral)
}

rcpp_nsk <- function(x, df, internal_knots, boundary_knots, trim = 0.0, complete_basis = FALSE, derivs = 0L, integral = FALSE) {
    .Call('_splines2_rcpp_nsk', PACKAGE = 'splines2', x, df, internal_knots, boundary_knots, trim, complete_basis, derivs, integral)
}

