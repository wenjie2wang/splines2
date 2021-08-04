// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rcpp_bSpline
Rcpp::NumericMatrix rcpp_bSpline(const arma::vec& x, const unsigned int df, const unsigned int degree, const arma::vec& internal_knots, const arma::vec& boundary_knots, const unsigned int derivs, const bool integral, const bool complete_basis);
RcppExport SEXP _splines2_rcpp_bSpline(SEXP xSEXP, SEXP dfSEXP, SEXP degreeSEXP, SEXP internal_knotsSEXP, SEXP boundary_knotsSEXP, SEXP derivsSEXP, SEXP integralSEXP, SEXP complete_basisSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type df(dfSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type degree(degreeSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type internal_knots(internal_knotsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type boundary_knots(boundary_knotsSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type derivs(derivsSEXP);
    Rcpp::traits::input_parameter< const bool >::type integral(integralSEXP);
    Rcpp::traits::input_parameter< const bool >::type complete_basis(complete_basisSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_bSpline(x, df, degree, internal_knots, boundary_knots, derivs, integral, complete_basis));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_bSpline_derivative
Rcpp::NumericMatrix rcpp_bSpline_derivative(const arma::vec& x, const unsigned int derivs, const unsigned int df, const unsigned int degree, const arma::vec& internal_knots, const arma::vec& boundary_knots, const bool complete_basis);
RcppExport SEXP _splines2_rcpp_bSpline_derivative(SEXP xSEXP, SEXP derivsSEXP, SEXP dfSEXP, SEXP degreeSEXP, SEXP internal_knotsSEXP, SEXP boundary_knotsSEXP, SEXP complete_basisSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type derivs(derivsSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type df(dfSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type degree(degreeSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type internal_knots(internal_knotsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type boundary_knots(boundary_knotsSEXP);
    Rcpp::traits::input_parameter< const bool >::type complete_basis(complete_basisSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_bSpline_derivative(x, derivs, df, degree, internal_knots, boundary_knots, complete_basis));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_bSpline_integral
Rcpp::NumericMatrix rcpp_bSpline_integral(const arma::vec& x, const unsigned int df, const unsigned int degree, const arma::vec& internal_knots, const arma::vec& boundary_knots, const bool complete_basis);
RcppExport SEXP _splines2_rcpp_bSpline_integral(SEXP xSEXP, SEXP dfSEXP, SEXP degreeSEXP, SEXP internal_knotsSEXP, SEXP boundary_knotsSEXP, SEXP complete_basisSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type df(dfSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type degree(degreeSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type internal_knots(internal_knotsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type boundary_knots(boundary_knotsSEXP);
    Rcpp::traits::input_parameter< const bool >::type complete_basis(complete_basisSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_bSpline_integral(x, df, degree, internal_knots, boundary_knots, complete_basis));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_bernsteinPoly
Rcpp::NumericMatrix rcpp_bernsteinPoly(const arma::vec& x, const unsigned int degree, const unsigned int derivs, const bool integral, const arma::vec& boundary_knots, const bool complete_basis);
RcppExport SEXP _splines2_rcpp_bernsteinPoly(SEXP xSEXP, SEXP degreeSEXP, SEXP derivsSEXP, SEXP integralSEXP, SEXP boundary_knotsSEXP, SEXP complete_basisSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type degree(degreeSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type derivs(derivsSEXP);
    Rcpp::traits::input_parameter< const bool >::type integral(integralSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type boundary_knots(boundary_knotsSEXP);
    Rcpp::traits::input_parameter< const bool >::type complete_basis(complete_basisSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_bernsteinPoly(x, degree, derivs, integral, boundary_knots, complete_basis));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_cSpline
Rcpp::NumericMatrix rcpp_cSpline(const arma::vec& x, const unsigned int df, const unsigned int degree, const arma::vec& internal_knots, const arma::vec& boundary_knots, const unsigned int derivs, const bool complete_basis);
RcppExport SEXP _splines2_rcpp_cSpline(SEXP xSEXP, SEXP dfSEXP, SEXP degreeSEXP, SEXP internal_knotsSEXP, SEXP boundary_knotsSEXP, SEXP derivsSEXP, SEXP complete_basisSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type df(dfSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type degree(degreeSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type internal_knots(internal_knotsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type boundary_knots(boundary_knotsSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type derivs(derivsSEXP);
    Rcpp::traits::input_parameter< const bool >::type complete_basis(complete_basisSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_cSpline(x, df, degree, internal_knots, boundary_knots, derivs, complete_basis));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_iSpline
Rcpp::NumericMatrix rcpp_iSpline(const arma::vec& x, const unsigned int df, const unsigned int degree, const arma::vec& internal_knots, const arma::vec& boundary_knots, const unsigned int derivs, const bool integral, const bool complete_basis);
RcppExport SEXP _splines2_rcpp_iSpline(SEXP xSEXP, SEXP dfSEXP, SEXP degreeSEXP, SEXP internal_knotsSEXP, SEXP boundary_knotsSEXP, SEXP derivsSEXP, SEXP integralSEXP, SEXP complete_basisSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type df(dfSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type degree(degreeSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type internal_knots(internal_knotsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type boundary_knots(boundary_knotsSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type derivs(derivsSEXP);
    Rcpp::traits::input_parameter< const bool >::type integral(integralSEXP);
    Rcpp::traits::input_parameter< const bool >::type complete_basis(complete_basisSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_iSpline(x, df, degree, internal_knots, boundary_knots, derivs, integral, complete_basis));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_mSpline
Rcpp::NumericMatrix rcpp_mSpline(const arma::vec& x, const unsigned int df, const unsigned int degree, const arma::vec& internal_knots, const arma::vec& boundary_knots, const unsigned int derivs, const bool integral, const bool complete_basis);
RcppExport SEXP _splines2_rcpp_mSpline(SEXP xSEXP, SEXP dfSEXP, SEXP degreeSEXP, SEXP internal_knotsSEXP, SEXP boundary_knotsSEXP, SEXP derivsSEXP, SEXP integralSEXP, SEXP complete_basisSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type df(dfSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type degree(degreeSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type internal_knots(internal_knotsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type boundary_knots(boundary_knotsSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type derivs(derivsSEXP);
    Rcpp::traits::input_parameter< const bool >::type integral(integralSEXP);
    Rcpp::traits::input_parameter< const bool >::type complete_basis(complete_basisSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_mSpline(x, df, degree, internal_knots, boundary_knots, derivs, integral, complete_basis));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_periodic_mSpline
Rcpp::NumericMatrix rcpp_periodic_mSpline(const arma::vec& x, const unsigned int df, const unsigned int degree, const arma::vec& internal_knots, const arma::vec& boundary_knots, const unsigned int derivs, const bool integral, const bool complete_basis);
RcppExport SEXP _splines2_rcpp_periodic_mSpline(SEXP xSEXP, SEXP dfSEXP, SEXP degreeSEXP, SEXP internal_knotsSEXP, SEXP boundary_knotsSEXP, SEXP derivsSEXP, SEXP integralSEXP, SEXP complete_basisSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type df(dfSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type degree(degreeSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type internal_knots(internal_knotsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type boundary_knots(boundary_knotsSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type derivs(derivsSEXP);
    Rcpp::traits::input_parameter< const bool >::type integral(integralSEXP);
    Rcpp::traits::input_parameter< const bool >::type complete_basis(complete_basisSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_periodic_mSpline(x, df, degree, internal_knots, boundary_knots, derivs, integral, complete_basis));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_naturalSpline
Rcpp::NumericMatrix rcpp_naturalSpline(const arma::vec& x, const unsigned int df, const arma::vec& internal_knots, const arma::vec& boundary_knots, const unsigned int derivs, const bool integral, const bool complete_basis);
RcppExport SEXP _splines2_rcpp_naturalSpline(SEXP xSEXP, SEXP dfSEXP, SEXP internal_knotsSEXP, SEXP boundary_knotsSEXP, SEXP derivsSEXP, SEXP integralSEXP, SEXP complete_basisSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type df(dfSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type internal_knots(internal_knotsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type boundary_knots(boundary_knotsSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type derivs(derivsSEXP);
    Rcpp::traits::input_parameter< const bool >::type integral(integralSEXP);
    Rcpp::traits::input_parameter< const bool >::type complete_basis(complete_basisSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_naturalSpline(x, df, internal_knots, boundary_knots, derivs, integral, complete_basis));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_splines2_rcpp_bSpline", (DL_FUNC) &_splines2_rcpp_bSpline, 8},
    {"_splines2_rcpp_bSpline_derivative", (DL_FUNC) &_splines2_rcpp_bSpline_derivative, 7},
    {"_splines2_rcpp_bSpline_integral", (DL_FUNC) &_splines2_rcpp_bSpline_integral, 6},
    {"_splines2_rcpp_bernsteinPoly", (DL_FUNC) &_splines2_rcpp_bernsteinPoly, 6},
    {"_splines2_rcpp_cSpline", (DL_FUNC) &_splines2_rcpp_cSpline, 7},
    {"_splines2_rcpp_iSpline", (DL_FUNC) &_splines2_rcpp_iSpline, 8},
    {"_splines2_rcpp_mSpline", (DL_FUNC) &_splines2_rcpp_mSpline, 8},
    {"_splines2_rcpp_periodic_mSpline", (DL_FUNC) &_splines2_rcpp_periodic_mSpline, 8},
    {"_splines2_rcpp_naturalSpline", (DL_FUNC) &_splines2_rcpp_naturalSpline, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_splines2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
