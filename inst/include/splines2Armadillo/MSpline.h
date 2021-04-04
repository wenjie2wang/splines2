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

#ifndef SPLINES2_MSPLINE_H
#define SPLINES2_MSPLINE_H

#include <stdexcept>

#include "common.h"
#include "utils.h"
#include "SplineBase.h"

namespace splines2 {

    // define a class for M-splines
    class MSpline : public SplineBase
    {
    protected:
        inline rmat get_basis_simple()
        {
            update_spline_df();
            update_x_index();
            update_knot_sequence();
            // define output matrix
            rmat b_mat {
                arma::zeros(x_.n_elem, spline_df_)
            };
            // generate basis of degree 0
            for (size_t i {0}; i < x_.n_elem; ++i) {
                unsigned int j { x_index_(i) };
                double denom { knot_sequence_(j + order_) -
                    knot_sequence_(j + degree_) };
                b_mat(i, j) = 1 / denom;
            }
            // main loop
            for (unsigned int k {1}; k <= degree_; ++k) {
                double dk { static_cast<double>(k) };
                double dk1 { (1.0 + 1.0 / dk) };
                const unsigned int k_offset { degree_ - k };
                // use the Cox-de Boor recursive formula
                for (size_t i {0}; i < x_.n_elem; ++i) {
                    double saved { 0 };
                    // for each x, at most "order" columns are not zero
                    // basis(x) is not zero from t_ii to t_{ii+k+1}
                    // where ii is index of x in terms of basis
                    for (size_t j {0}; j < k; ++j) {
                        size_t j_index { x_index_(i) + j };
                        size_t i1 { j_index + k_offset };
                        size_t i2 { j_index + order_ };
                        double den {
                            knot_sequence_(i2) - knot_sequence_(i1)
                        };
                        double term { dk1 * b_mat(i, j_index) };
                        b_mat(i, j_index) = saved +
                            (knot_sequence_(i2) - x_(i)) * term / den;
                        double den2 {
                            knot_sequence_(i2 + 1) -
                            knot_sequence_(i1 + 1)
                        };
                        saved = (x_(i) - knot_sequence_(i1 + 1)) *
                            term / den2;
                    }
                    b_mat(i, x_index_(i) + k) = saved;
                }
            }
            return b_mat;
        }

        inline rmat get_basis_extended()
        {
            MSpline msp_obj {
                x_, surrogate_internal_knots_, degree_,
                surrogate_boundary_knots_
            };
            rmat out { msp_obj.get_basis_simple() };
            // remove first and last #degree basis functions
            return out.cols(degree_, out.n_cols - order_);
        }

        inline rmat get_derivative_simple(
            const unsigned int derivs = 1
            )
        {
            MSpline msp_obj { this };
            msp_obj.set_degree(degree_ - derivs);
            // get basis matrix for (degree - derivs)
            rmat d_mat { msp_obj.get_basis_simple() };
            // make sure knot sequence and x index are latest
            update_knot_sequence();
            update_x_index();
            // add zero columns
            update_spline_df();
            // if (spline_df_ <= d_mat.n_cols) {
            //     throw std::range_error("FIXME: get_derivative_simple()");
            // }
            d_mat = add_zero_cols(d_mat, spline_df_ - d_mat.n_cols);
            // main loop
            for (unsigned int k {1}; k <= derivs; ++k) {
                const unsigned int k_offset { derivs - k };
                const size_t numer { degree_ - k_offset };
                for (size_t i {0}; i < x_.n_elem; ++i) {
                    double saved { 0 };
                    for (size_t j {0}; j < numer; ++j) {
                        size_t j_index { x_index_(i) + j };
                        size_t i1 { j_index + k_offset };
                        size_t i2 { j_index + order_ };
                        double den {
                            knot_sequence_(i2) - knot_sequence_(i1)
                        };
                        double term { (numer + 1) * d_mat(i, j_index) };
                        d_mat(i, j_index) = saved - term / den;
                        double den2 {
                            knot_sequence_(i2 + 1) - knot_sequence_(i1 + 1)
                        };
                        saved = term / den2;
                    }
                    d_mat(i, x_index_(i) + numer) = saved;
                }
            }
            return d_mat;
        }

        inline rmat get_derivative_extended(
            const unsigned int derivs = 1
            )
        {
            MSpline msp_obj {
                x_, surrogate_internal_knots_, degree_,
                surrogate_boundary_knots_
            };
            rmat out { msp_obj.get_derivative_simple(derivs) };
            // remove first and last #degree basis functions
            return out.cols(degree_, out.n_cols - order_);
        }

        inline rmat get_integral_simple()
        {
            // create a copy of this object
            MSpline msp_obj { this };
            // get basis matrix for (degree - derivs)
            msp_obj.set_degree(degree_ + 1);
            rmat i_mat { msp_obj.basis(false) };
            rvec knot_sequence_ord { msp_obj.get_knot_sequence() };
            // throw warning if any x is less than left-most boundary
            // if (arma::any(x_ < msp_obj.knot_sequence_(0))) {
            //     Rcpp::Rcout << "Warning: Found x < the leftmost knot, "
            //                 << msp_obj.knot_sequence_(0)
            //                 << ". "
            //                 << "The basis integrals were not well-defined."
            //                 << std::endl;
            // }
            // make sure x index is latest
            update_x_index();
            // compute t_{(i+1)+(k+1)+1} - t_{i+1} of s_{k}
            // which equals t_{(i+1)+(k+1)+2} - t_{(i+1)+1} of s_{k+1}
            arma::rowvec numer1 { arma::zeros<arma::rowvec>(i_mat.n_cols) };
            for (size_t j { 0 }; j < numer1.n_elem; ++j) {
                numer1(j) = knot_sequence_ord(j + order_ + 2) -
                    knot_sequence_ord(j + 1);
            }
            // for each row of i_mat
            for (size_t i {0}; i < x_.n_elem; ++i) {
                size_t k1 { x_index_(i) }, k2 { k1 + degree_ };
                arma::rowvec numer2 { i_mat(i, arma::span(k1, k2)) };
                arma::rowvec numer { numer1.cols(k1, k2) % numer2 };
                numer = rev_cum_sum(numer);
                for (size_t j {0}; j < i_mat.n_cols; ++j) {
                    if (j > k2) {
                        i_mat(i, j) = 0.0;
                    } else if (j >= k1) {
                        i_mat(i, j) = numer(j - k1) / (order_ + 1);
                    } else {
                        i_mat(i, j) = 1.0;
                    }
                }
            }
            return i_mat;
        }

        inline rmat get_integral_extended()
        {
            MSpline msp_obj {
                x_, surrogate_internal_knots_, degree_,
                surrogate_boundary_knots_
            };
            rmat out { msp_obj.get_integral_simple() };
            // remove first and last #degree basis functions
            return out.cols(degree_, out.n_cols - order_);
        }

    public:
        // inherits constructors
        using SplineBase::SplineBase;

        // function members

        //! Compute M-spline basis
        //!
        //! @param complete_basis A `bool` value indicating whether to return a
        //! complete spline basis
        //!
        //! @return arma::mat
        inline rmat basis(const bool complete_basis = true) override
        {
            rmat b_mat;
            if (is_extended_knot_sequence_) {
                b_mat = get_basis_extended();
            } else {
                b_mat = get_basis_simple();
            }
            // about to return
            if (complete_basis) {
                return b_mat;
            }
            // else
            return mat_wo_col1(b_mat);
        }

        // derivatives of M-splines
        inline rmat derivative(
            const unsigned int derivs = 1,
            const bool complete_basis = true
            ) override
        {
            if (derivs == 0) {
                throw std::range_error(
                    "'derivs' has to be a positive integer.");
            }
            // early exit if derivs is large enough
            update_spline_df();
            if (degree_ < derivs) {
                if (complete_basis) {
                    return arma::zeros(x_.n_elem, spline_df_);
                }
                if (spline_df_ == 1) {
                    throw std::range_error("No column left in the matrix.");
                }
                return arma::zeros(x_.n_elem, spline_df_ - 1);
            }
            rmat d_mat;
            if (is_extended_knot_sequence_) {
                d_mat = get_derivative_extended(derivs);
            } else {
                d_mat = get_derivative_simple(derivs);
            }
            // remove the first column if needed
            if (complete_basis) {
                return d_mat;
            }
            // else
            return mat_wo_col1(d_mat);
        }

        // integral of M-splines (I-splines)
        inline rmat integral(const bool complete_basis = true) override
        {
            rmat i_mat;
            if (is_extended_knot_sequence_) {
                i_mat = get_integral_extended();
            } else {
                i_mat = get_integral_simple();
            }
            // remove the first column if needed
            if (complete_basis) {
                return i_mat;
            }
            return mat_wo_col1(i_mat);
        }

    };
}  // splines2


#endif /* SPLINES2_MSPLINE_H */
