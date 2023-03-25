// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(splines2)]]
#include <splines2Armadillo.h>

// helper function
inline Rcpp::List return_list(splines2::NaturalSpline obj)
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
Rcpp::List rcpp_nspline00(const arma::vec& x,
                          const arma::vec& iknots,
                          const arma::vec& bknots)
{
    splines2::NaturalSpline obj;
    // setter methods
    obj.set_x(x)->
        set_internal_knots(iknots)->
        set_boundary_knots(bknots);
    return return_list(obj);
}
// [[Rcpp::export]]
Rcpp::List rcpp_nspline01(const arma::vec& x,
                          const arma::vec& iknots,
                          const arma::vec& bknots)
{
    splines2::NaturalSpline obj;
    // setter methods
    obj.set_x(x)->
        set_boundary_knots(bknots)->
        set_internal_knots(iknots);
    return return_list(obj);
}
// [[Rcpp::export]]
Rcpp::List rcpp_nspline02(const arma::vec& x,
                          const arma::vec& iknots,
                          const arma::vec& bknots)
{
    splines2::NaturalSpline obj;
    // setter methods
    obj.set_internal_knots(iknots)->
        set_x(x)->
        set_boundary_knots(bknots);
    return return_list(obj);
}
// [[Rcpp::export]]
Rcpp::List rcpp_nspline03(const arma::vec& x,
                          const arma::vec& iknots,
                          const arma::vec& bknots)
{
    splines2::NaturalSpline obj;
    // setter methods
    obj.set_boundary_knots(bknots)->
        set_x(x)->
        set_internal_knots(iknots);
    return return_list(obj);
}
// [[Rcpp::export]]
Rcpp::List rcpp_nspline04(const arma::vec& x,
                          const arma::vec& iknots,
                          const arma::vec& bknots)
{
    splines2::NaturalSpline obj;
    // setter methods
    obj.set_internal_knots(iknots)->
        set_boundary_knots(bknots)->
        set_x(x);
    return return_list(obj);
}
// [[Rcpp::export]]
Rcpp::List rcpp_nspline05(const arma::vec& x,
                          const arma::vec& iknots,
                          const arma::vec& bknots)
{
    splines2::NaturalSpline obj;
    // setter methods
    obj.set_boundary_knots(bknots)->
        set_internal_knots(iknots)->
        set_x(x);
    return return_list(obj);
}

// non-default constructor 1
// [[Rcpp::export]]
Rcpp::List rcpp_nspline1(const arma::vec& x,
                         const arma::vec& iknots,
                         const arma::vec& bknots)
{
    splines2::NaturalSpline obj { x, iknots, bknots };
    return return_list(obj);
}

// non-default constructor 2
// [[Rcpp::export]]
Rcpp::List rcpp_nspline2(const arma::vec& x,
                         const unsigned int df,
                         const arma::vec& bknots)
{
    splines2::NaturalSpline obj { x, df, bknots };
    return return_list(obj);
}

// non-default constructor 4
// [[Rcpp::export]]
Rcpp::List rcpp_nspline4(const arma::vec& x,
                         const arma::vec& iknots,
                         const arma::vec& bknots)
{
    splines2::NaturalSpline obj0;
    obj0.set_x(x)->
        set_internal_knots(iknots)->
        set_boundary_knots(bknots);
    splines2::NaturalSpline obj { &obj0 };
    return return_list(obj);
}

// conversion from BernsteinPoly
// [[Rcpp::export]]
Rcpp::List rcpp_nspline5(const arma::vec& x,
                         const arma::vec& iknots,
                         const unsigned int degree,
                         const arma::vec& bknots)
{
    splines2::BernsteinPoly obj0 { x, degree, bknots };
    splines2::NaturalSpline obj { static_cast<splines2::NaturalSpline>(obj0) };
    obj.set_internal_knots(iknots);
    return return_list(obj);
}

// conversion from PeriodicMSpline
// [[Rcpp::export]]
Rcpp::List rcpp_nspline6(const arma::vec& x,
                         const arma::vec& iknots,
                         const unsigned int degree,
                         const arma::vec& bknots)
{
    splines2::PeriodicMSpline obj0 { x, iknots, degree, bknots };
    splines2::NaturalSpline obj { static_cast<splines2::NaturalSpline>(obj0) };
    return return_list(obj);
}
