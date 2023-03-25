// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(splines2)]]
#include <splines2Armadillo.h>

// helper function
inline Rcpp::List return_list(splines2::BernsteinPoly obj)
{
    return Rcpp::List::create(
        Rcpp::Named("basis") = obj.basis(),
        Rcpp::Named("integral") = obj.integral(),
        Rcpp::Named("d1") = obj.derivative(),
        Rcpp::Named("d2") = obj.derivative(2),
        Rcpp::Named("d3") = obj.derivative(3),
        Rcpp::Named("degree") = obj.get_degree(),
        Rcpp::Named("boundary_knots") = splines2::arma2rvec(
            obj.get_boundary_knots())
        );
}


// default constructor and setter methods
// [[Rcpp::export]]
Rcpp::List rcpp_bp00(const arma::vec& x,
                     const unsigned int degree,
                     const arma::vec& bknots)
{
    splines2::BernsteinPoly obj;
    // setter methods
    obj.set_x(x)->
        set_boundary_knots(bknots)->
        set_degree(degree);
    return return_list(obj);
}
// [[Rcpp::export]]
Rcpp::List rcpp_bp01(const arma::vec& x,
                     const unsigned int degree,
                     const arma::vec& bknots)
{
    splines2::BernsteinPoly obj;
    // setter methods
    obj.set_x(x)->
        set_degree(degree)->
        set_boundary_knots(bknots);
    return return_list(obj);
}
// [[Rcpp::export]]
Rcpp::List rcpp_bp02(const arma::vec& x,
                     const unsigned int degree,
                     const arma::vec& bknots)
{
    splines2::BernsteinPoly obj;
    // setter methods
    obj.set_degree(degree)->
        set_x(x)->
        set_boundary_knots(bknots);
    return return_list(obj);
}
// [[Rcpp::export]]
Rcpp::List rcpp_bp03(const arma::vec& x,
                     const unsigned int degree,
                     const arma::vec& bknots)
{
    splines2::BernsteinPoly obj;
    // setter methods
    obj.set_degree(degree)->
        set_boundary_knots(bknots)->
        set_x(x);
    return return_list(obj);
}
// [[Rcpp::export]]
Rcpp::List rcpp_bp04(const arma::vec& x,
                     const unsigned int degree,
                     const arma::vec& bknots)
{
    splines2::BernsteinPoly obj;
    // setter methods
    obj.set_boundary_knots(bknots)->
        set_degree(degree)->
        set_x(x);
    return return_list(obj);
}
// [[Rcpp::export]]
Rcpp::List rcpp_bp05(const arma::vec& x,
                     const unsigned int degree,
                     const arma::vec& bknots)
{
    splines2::BernsteinPoly obj;
    // setter methods
    obj.set_boundary_knots(bknots)->
        set_x(x)->
        set_degree(degree);
    return return_list(obj);
}

// non-default constructor 1
// [[Rcpp::export]]
Rcpp::List rcpp_bp1(const arma::vec& x,
                    const unsigned int degree,
                    const arma::vec& bknots)
{
    splines2::BernsteinPoly obj { x, degree, bknots };
    return return_list(obj);
}

// non-default constructor 4
// [[Rcpp::export]]
Rcpp::List rcpp_bp4(const arma::vec& x,
                    const arma::vec& iknots,
                    const unsigned int degree,
                    const arma::vec& bknots)
{
    splines2::MSpline obj0;
    obj0.set_x(x)->
        set_degree(degree)->
        set_internal_knots(iknots)->
        set_boundary_knots(bknots);
    splines2::BernsteinPoly obj { &obj0 };
    return return_list(obj);
}

// conversion from PeriodicMSpline
// [[Rcpp::export]]
Rcpp::List rcpp_bp6(const arma::vec& x,
                    const arma::vec& iknots,
                    const unsigned int degree,
                    const arma::vec& bknots)
{
    splines2::PeriodicMSpline obj0 { x, iknots, degree, bknots };
    splines2::BernsteinPoly obj {
        static_cast<splines2::BernsteinPoly>(obj0)
    };
    return return_list(obj);
}
