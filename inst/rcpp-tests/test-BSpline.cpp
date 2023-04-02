// [[Rcpp::depends(RcppArmadillo)]]
#include "../include/splines2Armadillo.h"

// helper function
inline Rcpp::List return_list(splines2::BSpline obj)
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
Rcpp::List rcpp_bspline00(const arma::vec& x,
                          const arma::vec& iknots,
                          const unsigned int degree,
                          const arma::vec& bknots)
{
    splines2::BSpline obj;
    // setter methods
    obj.set_x(x)->
        set_internal_knots(iknots)->
        set_boundary_knots(bknots)->
        set_degree(degree);
    return return_list(obj);
}
// [[Rcpp::export]]
Rcpp::List rcpp_bspline01(const arma::vec& x,
                          const arma::vec& iknots,
                          const unsigned int degree,
                          const arma::vec& bknots)
{
    splines2::BSpline obj;
    // setter methods
    obj.set_x(x)->
        set_boundary_knots(bknots)->
        set_internal_knots(iknots)->
        set_degree(degree);
    return return_list(obj);
}
// [[Rcpp::export]]
Rcpp::List rcpp_bspline02(const arma::vec& x,
                          const arma::vec& iknots,
                          const unsigned int degree,
                          const arma::vec& bknots)
{
    splines2::BSpline obj;
    // setter methods
    obj.set_internal_knots(iknots)->
        set_x(x)->
        set_boundary_knots(bknots)->
        set_degree(degree);
    return return_list(obj);
}
// [[Rcpp::export]]
Rcpp::List rcpp_bspline03(const arma::vec& x,
                          const arma::vec& iknots,
                          const unsigned int degree,
                          const arma::vec& bknots)
{
    splines2::BSpline obj;
    // setter methods
    obj.set_degree(degree)->
        set_boundary_knots(bknots)->
        set_x(x)->
        set_internal_knots(iknots);
    return return_list(obj);
}
// [[Rcpp::export]]
Rcpp::List rcpp_bspline04(const arma::vec& x,
                          const arma::vec& iknots,
                          const unsigned int degree,
                          const arma::vec& bknots)
{
    splines2::BSpline obj;
    // setter methods
    obj.set_degree(degree)->
        set_internal_knots(iknots)->
        set_boundary_knots(bknots)->
        set_x(x);
    return return_list(obj);
}
// [[Rcpp::export]]
Rcpp::List rcpp_bspline05(const arma::vec& x,
                          const arma::vec& iknots,
                          const unsigned int degree,
                          const arma::vec& bknots)
{
    splines2::BSpline obj;
    // setter methods
    obj.set_degree(degree)->
        set_boundary_knots(bknots)->
        set_internal_knots(iknots)->
        set_x(x);
    return return_list(obj);
}

// non-default constructor 1
// [[Rcpp::export]]
Rcpp::List rcpp_bspline1(const arma::vec& x,
                         const arma::vec& iknots,
                         const unsigned int degree,
                         const arma::vec& bknots)
{
    splines2::BSpline obj { x, iknots, degree, bknots };
    return return_list(obj);
}

// non-default constructor 2
// [[Rcpp::export]]
Rcpp::List rcpp_bspline2(const arma::vec& x,
                         const unsigned int df,
                         const unsigned int degree,
                         const arma::vec& bknots)
{
    splines2::BSpline obj { x, df, degree, bknots };
    return return_list(obj);
}

// non-default constructor 3
// [[Rcpp::export]]
Rcpp::List rcpp_bspline3(const arma::vec& x,
                         const unsigned int degree,
                         const arma::vec& knot_seq)
{
    splines2::BSpline obj { x, degree, knot_seq };
    return return_list(obj);
}

// non-default constructor 4
// [[Rcpp::export]]
Rcpp::List rcpp_bspline4(const arma::vec& x,
                         const arma::vec& iknots,
                         const unsigned int degree,
                         const arma::vec& bknots)
{
    splines2::MSpline obj0;
    obj0.set_x(x)->
        set_internal_knots(iknots)->
        set_degree(degree)->
        set_boundary_knots(bknots);
    splines2::BSpline obj { &obj0 };
    return return_list(obj);
}

// conversion from BernsteinPoly
// [[Rcpp::export]]
Rcpp::List rcpp_bspline5(const arma::vec& x,
                         const arma::vec& iknots,
                         const unsigned int degree,
                         const arma::vec& bknots)
{
    splines2::BernsteinPoly obj0 { x, degree, bknots };
    splines2::BSpline obj { static_cast<splines2::BSpline>(obj0) };
    obj.set_internal_knots(iknots);
    return return_list(obj);
}

// conversion from PeriodicMSpline
// [[Rcpp::export]]
Rcpp::List rcpp_bspline6(const arma::vec& x,
                         const arma::vec& iknots,
                         const unsigned int degree,
                         const arma::vec& bknots)
{
    splines2::PeriodicMSpline obj0 { x, iknots, degree, bknots };
    splines2::BSpline obj { static_cast<splines2::BSpline>(obj0) };
    return return_list(obj);
}
