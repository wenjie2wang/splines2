// [[Rcpp::depends(RcppArmadillo)]]
#include "../include/splines2Armadillo.h"

// helper function
inline Rcpp::List return_list(splines2::PeriodicMSpline obj)
{
    return Rcpp::List::create(
        Rcpp::Named("basis") = obj.basis(),
        Rcpp::Named("integral") = obj.integral(),
        Rcpp::Named("d1") = obj.derivative(),
        Rcpp::Named("d2") = obj.derivative(2),
        Rcpp::Named("d3") = obj.derivative(3),
        Rcpp::Named("degree") = obj.get_degree(),
        Rcpp::Named("internal_knots") = splines2::arma2rvec(
            obj.get_internal_knots()),
        Rcpp::Named("boundary_knots") = splines2::arma2rvec(
            obj.get_boundary_knots())
        );
}


// default constructor and setter methods
// [[Rcpp::export]]
Rcpp::List rcpp_pmspline00(const arma::vec& x,
                           const arma::vec& iknots,
                           const unsigned int degree,
                           const arma::vec& bknots)
{
    splines2::PeriodicMSpline obj;
    obj.set_x(x)->
        set_internal_knots(iknots)->
        set_boundary_knots(bknots)->
        set_degree(degree);
    return return_list(obj);
}
// [[Rcpp::export]]
Rcpp::List rcpp_pmspline01(const arma::vec& x,
                           const arma::vec& iknots,
                           const unsigned int degree,
                           const arma::vec& bknots)
{
    splines2::PeriodicMSpline obj;
    obj.set_x(x)->
        set_boundary_knots(bknots)->
        set_internal_knots(iknots)->
        set_degree(degree);
    return return_list(obj);
}
// [[Rcpp::export]]
Rcpp::List rcpp_pmspline02(const arma::vec& x,
                           const arma::vec& iknots,
                           const unsigned int degree,
                           const arma::vec& bknots)
{
    splines2::PeriodicMSpline obj;
    obj.set_internal_knots(iknots)->
        set_x(x)->
        set_boundary_knots(bknots)->
        set_degree(degree);
    return return_list(obj);
}
// [[Rcpp::export]]
Rcpp::List rcpp_pmspline03(const arma::vec& x,
                           const arma::vec& iknots,
                           const unsigned int degree,
                           const arma::vec& bknots)
{
    splines2::PeriodicMSpline obj;
    obj.set_degree(degree)->
        set_boundary_knots(bknots)->
        set_x(x)->
        set_internal_knots(iknots);
    return return_list(obj);
}
// [[Rcpp::export]]
Rcpp::List rcpp_pmspline04(const arma::vec& x,
                           const arma::vec& iknots,
                           const unsigned int degree,
                           const arma::vec& bknots)
{
    splines2::PeriodicMSpline obj;
    obj.set_degree(degree)->
        set_internal_knots(iknots)->
        set_boundary_knots(bknots)->
        set_x(x);
    return return_list(obj);
}
// [[Rcpp::export]]
Rcpp::List rcpp_pmspline05(const arma::vec& x,
                           const arma::vec& iknots,
                           const unsigned int degree,
                           const arma::vec& bknots)
{
    splines2::PeriodicMSpline obj;
    obj.set_degree(degree)->
        set_boundary_knots(bknots)->
        set_internal_knots(iknots)->
        set_x(x);
    return return_list(obj);
}
// [[Rcpp::export]]
Rcpp::List rcpp_pmspline06(const arma::vec& x,
                           const unsigned int degree,
                           const arma::vec& knot_seq)
{
    splines2::PeriodicMSpline obj;
    obj.set_degree(degree)->
        set_knot_sequence(knot_seq)->
        set_x(x);
    return return_list(obj);
}
// [[Rcpp::export]]
Rcpp::List rcpp_pmspline07(const arma::vec& x,
                           const unsigned int degree,
                           const arma::vec& knot_seq)
{
    splines2::PeriodicMSpline obj;
    obj.set_knot_sequence(knot_seq)->
        set_degree(degree)->
        set_x(x);
    return return_list(obj);
}
// [[Rcpp::export]]
Rcpp::List rcpp_pmspline08(const arma::vec& x,
                           const unsigned int degree,
                           const arma::vec& knot_seq)
{
    splines2::PeriodicMSpline obj;
    obj.set_x(x)->
        set_degree(degree)->
        set_knot_sequence(knot_seq);
    return return_list(obj);
}

// non-default constructor 1
// [[Rcpp::export]]
Rcpp::List rcpp_pmspline1(const arma::vec& x,
                          const arma::vec& iknots,
                          const unsigned int degree,
                          const arma::vec& bknots)
{
    splines2::PeriodicMSpline obj { x, iknots, degree, bknots };
    return return_list(obj);
}

// non-default constructor 2
// [[Rcpp::export]]
Rcpp::List rcpp_pmspline2(const arma::vec& x,
                          const unsigned int df,
                          const unsigned int degree,
                          const arma::vec& bknots)
{
    splines2::PeriodicMSpline obj { x, df, degree, bknots };
    return return_list(obj);
}

// non-default constructor 3
// [[Rcpp::export]]
Rcpp::List rcpp_pmspline3(const arma::vec& x,
                          const unsigned int degree,
                          const arma::vec& knot_seq)
{
    splines2::PeriodicMSpline obj { x, degree, knot_seq };
    return return_list(obj);
}

// non-default constructor 4
// [[Rcpp::export]]
Rcpp::List rcpp_pmspline4(const arma::vec& x,
                          const arma::vec& iknots,
                          const unsigned int degree,
                          const arma::vec& bknots)
{
    splines2::MSpline obj0;
    obj0.set_x(x)->
        set_internal_knots(iknots)->
        set_degree(degree)->
        set_boundary_knots(bknots);
    splines2::PeriodicMSpline obj { &obj0 };
    return return_list(obj);
}

// conversion from BernsteinPoly
// [[Rcpp::export]]
Rcpp::List rcpp_pmspline5(const arma::vec& x,
                          const arma::vec& iknots,
                          const unsigned int degree,
                          const arma::vec& bknots)
{
    splines2::BernsteinPoly obj0 { x, degree, bknots };
    splines2::PeriodicMSpline obj {
        static_cast<splines2::PeriodicMSpline>(obj0)
    };
    obj.set_internal_knots(iknots);
    return return_list(obj);
}
