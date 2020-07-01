#include <RcppArmadillo.h>
// [[Rcpp::plugins(cpp11)]]

#include <splines2Armadillo.h>

// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_bSpline(const arma::vec& x,
                                 const unsigned int df,
                                 const unsigned int degree,
                                 const arma::vec& internal_knots,
                                 const arma::vec& boundary_knots,
                                 const bool complete_basis = true)
{
    const unsigned int wo_intercept {
        static_cast<unsigned int>(! complete_basis)
    };
    splines2::BSpline bs_obj;
    // if df > 0 and knots are not specified
    // auto set internal knots based on df
    if (df > 0 && internal_knots.n_elem == 0) {
        // compute actual spline degree of freedom
        unsigned int spline_df { df + wo_intercept };
        bs_obj = splines2::BSpline(x, spline_df, degree, boundary_knots);
    } else {
        // else ignore df
        bs_obj = splines2::BSpline(x, internal_knots, degree, boundary_knots);
    }
    Rcpp::NumericMatrix out {
        splines2::arma2rmat(bs_obj.basis(complete_basis))
            };
    // add attributes
    out.attr("dimnames") = Rcpp::List::create(
        R_NilValue, splines2::char_seq_len(out.ncol())
        );
    out.attr("x") = splines2::arma2rvec(x);
    out.attr("degree") = bs_obj.get_degree();
    out.attr("knots") = splines2::arma2rvec(bs_obj.get_internal_knots());
    out.attr("Boundary.knots") =
        splines2::arma2rvec(bs_obj.get_boundary_knots());
    out.attr("intercept") = complete_basis;
    return out;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_dbs(const arma::vec& x,
                             const unsigned int derivs,
                             const unsigned int df,
                             const unsigned int degree,
                             const arma::vec& internal_knots,
                             const arma::vec& boundary_knots,
                             const bool complete_basis = true)
{
    const unsigned int wo_intercept {
        static_cast<unsigned int>(! complete_basis)
    };
    splines2::BSpline bs_obj;
    // if df > 0 and knots are not specified
    // auto set internal knots based on df
    if (df > 0 && internal_knots.n_elem == 0) {
        // compute actual spline degree of freedom
        unsigned int spline_df { df + wo_intercept };
        bs_obj = splines2::BSpline(x, spline_df, degree, boundary_knots);
    } else {
        // else ignore df
        bs_obj = splines2::BSpline(x, internal_knots, degree, boundary_knots);
    }
    Rcpp::NumericMatrix out {
        splines2::arma2rmat(bs_obj.derivative(complete_basis))
            };
    // add attributes
    out.attr("dimnames") = Rcpp::List::create(
        R_NilValue, splines2::char_seq_len(out.ncol())
        );
    out.attr("x") = splines2::arma2rvec(x);
    out.attr("derivs") = derivs;
    out.attr("degree") = bs_obj.get_degree();
    out.attr("knots") = splines2::arma2rvec(bs_obj.get_internal_knots());
    out.attr("Boundary.knots") =
        splines2::arma2rvec(bs_obj.get_boundary_knots());
    out.attr("intercept") = complete_basis;
    return out;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_ibs(const arma::vec& x,
                             const unsigned int df,
                             const unsigned int degree,
                             const arma::vec& internal_knots,
                             const arma::vec& boundary_knots,
                             const bool complete_basis = true)
{
    const unsigned int wo_intercept {
        static_cast<unsigned int>(! complete_basis)
    };
    splines2::BSpline bs_obj;
    // if df > 0 and knots are not specified
    // auto set internal knots based on df
    if (df > 0 && internal_knots.n_elem == 0) {
        // compute actual spline degree of freedom
        unsigned int spline_df { df + wo_intercept };
        bs_obj = splines2::BSpline(x, spline_df, degree, boundary_knots);
    } else {
        // else ignore df
        bs_obj = splines2::BSpline(x, internal_knots, degree, boundary_knots);
    }
    Rcpp::NumericMatrix out {
        splines2::arma2rmat(bs_obj.integral(complete_basis))
            };
    // add attributes
    out.attr("dimnames") = Rcpp::List::create(
        R_NilValue, splines2::char_seq_len(out.ncol())
        );
    out.attr("x") = splines2::arma2rvec(x);
    out.attr("degree") = bs_obj.get_degree();
    out.attr("knots") = splines2::arma2rvec(bs_obj.get_internal_knots());
    out.attr("Boundary.knots") =
        splines2::arma2rvec(bs_obj.get_boundary_knots());
    out.attr("intercept") = complete_basis;
    return out;
}
