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

#ifndef SPLINES2_BSPLINE_H
#define SPLINES2_BSPLINE_H

#include <stdexcept>

#include "common.h"
#include "utils.h"
#include "SplineBase.h"

namespace splines2 {

    // define a class for B-splines
    class BSpline : public SplineBase
    {
    protected:
        inline virtual rmat get_basis_simple()
        {
            update_spline_df();
            update_x_index();
            // define output matrix
            rmat b_mat {
                arma::zeros(x_.n_elem, spline_df_)
            };
            // generate basis of degree 0
            for (size_t i {0}; i < x_.n_elem; ++i) {
                b_mat(i, x_index_(i)) = 1;
            }
            // only need knot sequence for degree > 0 unless extended
            if (degree_ > 0) {
                update_knot_sequence();
            }
            if (has_internal_multiplicity_) {
                for (unsigned int k {1}; k <= degree_; ++k) {
                    const unsigned int k_offset { degree_ - k };
                    // use the Cox-de Boor recursive formula
                    for (size_t i {0}; i < x_.n_elem; ++i) {
                        double saved { 0.0 };
                        for (size_t j {0}; j < k; ++j) {
                            size_t j_index { x_index_(i) + j };
                            size_t i1 { j_index + k_offset + 1 };
                            size_t i2 { j_index + order_ };
                            double den {
                                knot_sequence_(i2) - knot_sequence_(i1)
                            };
                            if (isAlmostEqual(den)) {
                                if (j != 0 || knot_sequence_(i2) - x_(i) != 0) {
                                    b_mat(i, j_index) = saved;
                                }
                                saved = 0.0;
                            } else {
                                double term { b_mat(i, j_index) / den };
                                b_mat(i, j_index) = saved +
                                    (knot_sequence_(i2) - x_(i)) * term;
                                saved = (x_(i) - knot_sequence_(i1)) * term;
                            }
                        }
                        b_mat(i, x_index_(i) + k) = saved;
                    }
                }
            } else {
                for (unsigned int k {1}; k <= degree_; ++k) {
                    const unsigned int k_offset { degree_ - k };
                    // use the Cox-de Boor recursive formula
                    for (size_t i {0}; i < x_.n_elem; ++i) {
                        double saved { 0.0 };
                        // for each x, at most "order" columns are not zero
                        // basis_j(x) is not zero from t_j to t_{j+k+1}
                        // where j is index of x in terms of basis
                        // knot sequence: t0, t1, t2, ...
                        for (size_t j {0}; j < k; ++j) {
                            size_t j_index { x_index_(i) + j };
                            size_t i1 { j_index + k_offset + 1 };
                            size_t i2 { j_index + order_ };
                            double den {
                                knot_sequence_(i2) - knot_sequence_(i1)
                            };
                            // no need to check for distinct internal knots
                            double term { b_mat(i, j_index) / den };
                            b_mat(i, j_index) = saved +
                                (knot_sequence_(i2) - x_(i)) * term;
                            saved = (x_(i) - knot_sequence_(i1)) * term;
                        }
                        b_mat(i, x_index_(i) + k) = saved;
                    }
                }
            }
            return b_mat;
        }

        inline virtual rmat get_basis_extended()
        {
            BSpline bsp_obj {
                x_, surrogate_internal_knots_, degree_,
                surrogate_boundary_knots_
            };
            rmat out { bsp_obj.basis(true) };
            // remove first and last #degree basis functions
            return out.cols(degree_, out.n_cols - order_);
        }

        inline virtual rmat get_derivative_simple(
            const unsigned int derivs = 1
            )
        {
            // create a copy of this object
            BSpline bsp_obj { this };
            bsp_obj.set_degree(degree_ - derivs);
            // get basis matrix for (degree - derivs)
            rmat d_mat { bsp_obj.basis(true) };
            // make sure knot sequence and x index are latest
            update_knot_sequence();
            update_x_index();
            // add zero columns
            d_mat = add_zero_cols(d_mat, spline_df_ - d_mat.n_cols);
            if (has_internal_multiplicity_) {
                for (unsigned int k {1}; k <= derivs; ++k) {
                    const unsigned int k_offset { derivs - k };
                    const size_t numer { degree_ - k_offset };
                    for (size_t i {0}; i < x_.n_elem; ++i) {
                        double saved { 0 };
                        for (size_t j {0}; j < numer; ++j) {
                            size_t j_index { x_index_(i) + j };
                            size_t i1 { j_index + k_offset + 1 };
                            size_t i2 { j_index + order_ };
                            double den {
                                knot_sequence_(i2) - knot_sequence_(i1)
                            };
                            if (isAlmostEqual(den)) {
                                if (j != 0 || knot_sequence_(i2) - x_(i) != 0) {
                                    d_mat(i, j_index) = saved;
                                }
                                saved = 0.0;
                            } else {
                                double term { numer * d_mat(i, j_index) / den };
                                d_mat(i, j_index) = saved - term;
                                saved = term;
                            }
                        }
                        d_mat(i, x_index_(i) + numer) = saved;
                    }
                }
            } else {
                for (unsigned int k {1}; k <= derivs; ++k) {
                    const unsigned int k_offset { derivs - k };
                    const size_t numer { degree_ - k_offset };
                    for (size_t i {0}; i < x_.n_elem; ++i) {
                        double saved { 0 };
                        for (size_t j {0}; j < numer; ++j) {
                            size_t j_index { x_index_(i) + j };
                            size_t i1 { j_index + k_offset + 1 };
                            size_t i2 { j_index + order_ };
                            double den {
                                knot_sequence_(i2) - knot_sequence_(i1)
                            };
                            double term { numer * d_mat(i, j_index) / den };
                            d_mat(i, j_index) = saved - term;
                            saved = term;
                        }
                        d_mat(i, x_index_(i) + numer) = saved;
                    }
                }
            }
            return d_mat;
        }

        inline virtual rmat get_derivative_extended(
            const unsigned int derivs = 1
            )
        {
            BSpline bsp_obj {
                x_, surrogate_internal_knots_, degree_,
                surrogate_boundary_knots_
            };
            rmat out { bsp_obj.derivative(derivs, true) };
            // remove first and last #degree basis functions
            return out.cols(degree_, out.n_cols - order_);
        }

        inline virtual rmat get_integral_simple()
        {
            // create a copy of this object
            BSpline bsp_obj { this };
            bsp_obj.set_degree(degree_ + 1);
            rmat i_mat { bsp_obj.basis(false) };
            rvec knot_sequence_ord { bsp_obj.get_knot_sequence() };
            // make sure x index are latest
            update_x_index();
            // compute t_{(i+1)+k+1} - t_{i+1} of s_{k}
            // which is t_{(i+1)+(k+1)+1} - t_{(i+1)+1} of s_{k+1}
            arma::rowvec numer1 { arma::zeros<arma::rowvec>(i_mat.n_cols) };
            for (size_t j { 0 }; j < numer1.n_elem; ++j) {
                numer1(j) = knot_sequence_ord(j + order_ + 1) -
                    knot_sequence_ord(j + 1);
            }
            // for each row of i_mat
            for (size_t i {0}; i < x_.n_elem; ++i) {
                size_t k1 { x_index_(i) }, k2 { k1 + degree_ };
                arma::rowvec numer2 { i_mat(i, arma::span(k1, k2)) };
                numer2 = rev_cum_sum(numer2);
                for (size_t j {0}; j < i_mat.n_cols; ++j) {
                    if (j > k2) {
                        i_mat(i, j) = 0;
                    } else if (j >= k1) {
                        i_mat(i, j) = numer2(j - k1) * numer1(j) / order_;
                    } else {
                        i_mat(i, j) = numer1(j) / order_;
                    }
                }
            }
            return i_mat;
        }

        inline virtual rmat get_integral_extended()
        {
            BSpline bsp_obj {
                x_, surrogate_internal_knots_, degree_,
                surrogate_boundary_knots_
            };
            rmat out { bsp_obj.integral(true) };
            // remove first and last #degree basis functions
            return out.cols(degree_, out.n_cols - order_);
        }

    public:
        // inherits constructors
        using SplineBase::SplineBase;

        // function members

        //! Compute B-spline basis functions
        //!
        //! @param complete_basis A `bool` value indicating whether to return a
        //! complete spline basis
        //!
        //! @return arma::mat
        inline virtual rmat basis(const bool complete_basis = true)
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

        // derivatives of B-splines
        inline virtual rmat derivative(
            const unsigned int derivs = 1,
            const bool complete_basis = true
            )
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

        // integral of B-splines
        inline virtual rmat integral(const bool complete_basis = true)
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

}

#endif
