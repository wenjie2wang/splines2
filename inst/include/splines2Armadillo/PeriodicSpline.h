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

#ifndef SPLINES2_PERIODICSPLINE_H
#define SPLINES2_PERIODICSPLINE_H

#include <stdexcept>

#include "common.h"
#include "utils.h"
#include "SplineBase.h"
#include "MSpline.h"
#include "BSpline.h"

namespace splines2 {

    // define a class for nonnegative periodic spline basis over [a, b]
    // with unit integral over [a, b]
    template <typename T_sp>
    class PeriodicSpline : public SplineBase
    {
    protected:
        double range_size_;     // b - a
        rvec x_in_range_;       // x in range
        rvec x_num_shift_;      // x_ = x_in_range + x_num_shift * range_size
        bool is_x_in_range_latest_ = false;

        // compute range size given boundary knots
        inline void update_range_size()
        {
            range_size_ = boundary_knots_(1) - boundary_knots_(0);
        }
        // compute spline df
        inline void update_spline_df() override
        {
            spline_df_ = internal_knots_.n_elem + 1;
        }

        // extend knot sequence for periodic splines
        // reference:
        // - Farin, G., & Hansford, D. (2000). The essentials of CAGD.
        // - Schumaker, L. (2007). Spline functions: Basic theory.
        // - splines::splineDesign() when outer.ok = TRUE
        // - Piegl, L., & Tiller, W. (1997). The NURBS book. Springer
        inline void extend_knot_sequence()
        {
            // number of internal knots must >= degree
            if (internal_knots_.n_elem + 1 < degree_) {
                throw std::range_error(
                    "The number of unique internal knots must be >= degree - 1."
                    );
            }
            // extend to distinct knot sequence first
            rvec res { arma::zeros(internal_knots_.n_elem + 2 * order_) };
            // place the internal knots
            for (size_t i {0}; i < internal_knots_.n_elem; ++i) {
                res(degree_ + 1 + i) = internal_knots_(i);
            }
            // place the boundary knots
            res(degree_) = boundary_knots_(0);
            res(res.n_elem - 1 - degree_) = boundary_knots_(1);
            // first and last "degree-size" knots

            size_t n_ { degree_ + internal_knots_.n_elem };
            for (size_t i {0}; i < degree_; ++i) {
                res(degree_ - i - 1) = res(degree_ - i) -
                    (res(n_ - i + 1) - res(n_ - i));
                res(n_ + i + 2) = res(n_ + i + 1) +
                    (res(degree_ + i + 1) - res(degree_ + i));
            }
            surrogate_boundary_knots_ = arma::zeros(2);
            surrogate_boundary_knots_(0) = res(0);
            surrogate_boundary_knots_(1) = res(res.n_elem - 1);
            surrogate_internal_knots_ = res.subvec(1, res.n_elem - 2);
            // add multiplicities to "boundary" knots
            knot_sequence_ = arma::zeros(res.n_elem + 2 * degree_);
            for (size_t i {0}; i < knot_sequence_.n_elem; ++i) {
                if (i < order_) {
                    knot_sequence_(i) = surrogate_boundary_knots_(0);
                } else if (i < knot_sequence_.n_elem - order_) {
                    knot_sequence_(i) = surrogate_internal_knots_(i - order_);
                } else {
                    knot_sequence_(i) = surrogate_boundary_knots_(1);
                }
            }
        }

        inline void update_knot_sequence() override
        {
            if (! is_knot_sequence_latest_ || knot_sequence_.n_elem == 0) {
                if (is_extended_knot_sequence_) {
                    set_extended_knot_sequence(knot_sequence_);
                } else {
                    // if it is a simple knot sequence previously,
                    // we will assume a simple knot sequence
                    set_simple_knot_sequence();
                }
            }
            stopifnot_simple_knot_sequence();
            extend_knot_sequence();
        }

        inline void set_x_in_range()
        {
            if (is_x_in_range_latest_) {
                return;
            }
            update_range_size();
            x_num_shift_ = arma::floor((x_ - boundary_knots_(0)) / range_size_);
            x_in_range_ = x_ - range_size_ * x_num_shift_;
        }

        inline rmat clamp_basis(const rmat& b_mat)
        {
            rmat out { b_mat.head_cols(degree_) + b_mat.tail_cols(degree_) };
            if (internal_knots_.n_elem + 1 > degree_) {
                rmat out1 { b_mat.cols(degree_, internal_knots_.n_elem) };
                out = arma::join_rows(out1, out);
            }
            return out;
        }

    public:
        PeriodicSpline() {}

        explicit PeriodicSpline(const SplineBase* pSplineBase) :
            SplineBase(pSplineBase)
        {
            // stopifnot_simple_knot_sequence();
            update_spline_df();
        }

        // explicit conversion
        template <typename T>
        explicit operator T() const {
            T obj;
            obj.set_x(x_)->
                set_degree(degree_)->
                set_internal_knots(internal_knots_)->
                set_boundary_knots(boundary_knots_);
            return obj;
        }

        // given boundary_knots for consistency with SplineBase
        PeriodicSpline<T_sp>(const rvec& x,
                             const rvec& internal_knots,
                             const unsigned int degree = 3,
                             const rvec& boundary_knots = rvec())
        {
            x_ = x;
            degree_ = degree;
            simplify_knots(internal_knots, boundary_knots);
            if (internal_knots_.n_elem + 1 < degree_) {
                throw std::range_error(
                    "The number of unique internal knots must be >= degree - 1."
                    );
            }
            order_ = degree_ + 1;
            update_spline_df();
        }

        PeriodicSpline<T_sp>(const rvec& x,
                             const unsigned int spline_df,
                             const unsigned int degree = 3,
                             const rvec& boundary_knots = rvec())
        {
            x_ = x;
            degree_ = degree;
            // spline_df = number(internal_knot) + 1 >= degree
            if (spline_df < degree) {
                throw std::range_error(
                    "The specified 'df' must be > 'degree'.");
            }
            spline_df_ = spline_df;
            order_ = degree_ + 1;
            // determine internal knots by spline_df and x
            unsigned int n_internal_knots { spline_df_ - 1 };
            rvec prob_vec { arma::linspace(0, 1, n_internal_knots + 2) };
            prob_vec = prob_vec.subvec(1, n_internal_knots);
            simplify_knots(rvec(), boundary_knots);
            // get quantiles of x in range
            set_x_in_range();
            rvec internal_knots { arma_quantile(x_in_range_, prob_vec) };
            simplify_knots(internal_knots);
        }

        // possible to specify knot sequence directly. But it must be a simple
        // knot sequence.
        PeriodicSpline<T_sp>(const rvec& x,
                             const unsigned int degree,
                             const rvec& knot_sequence)
        {
            x_ = x;
            degree_ = degree;
            order_ = degree_ + 1;
            set_extended_knot_sequence(knot_sequence);
            stopifnot_simple_knot_sequence();
        }
        inline SplineBase* set_knot_sequence(
            const rvec& knot_sequence
            ) override
        {
            set_extended_knot_sequence(knot_sequence);
            return this;
        }

        //! Compute periodic spline basis based on M-spline basis
        //!
        //! @param complete_basis A `bool` value indicating whether to return a
        //! complete spline basis
        //!
        //! @return arma::mat
        inline rmat basis(const bool complete_basis = true) override
        {
            update_knot_sequence();
            set_x_in_range();
            // create a Spline object for the extended knot sequence
            T_sp bs_obj { x_in_range_, degree_, knot_sequence_ };
            rmat b_mat { bs_obj.basis(true) };
            // remove first and last #degree basis functions
            b_mat = b_mat.cols(degree_, b_mat.n_cols - order_);
            b_mat = clamp_basis(b_mat);
            // post-processing
            if (complete_basis) {
                return b_mat;
            }
            // else
            return mat_wo_col1(b_mat);
        }

        inline rmat derivative(
            const unsigned int derivs = 1,
            const bool complete_basis = true
            ) override
        {
            update_knot_sequence();
            if (derivs == 0) {
                throw std::range_error(
                    "'derivs' has to be a positive integer.");
            }
            // early exit if derivs is large enough
            if (degree_ < derivs) {
                unsigned int old_df { spline_df_ };
                if (complete_basis) {
                    return arma::zeros(x_.n_elem, old_df);
                }
                if (old_df == 1) {
                    throw std::range_error("No column left in the matrix.");
                }
                return arma::zeros(x_.n_elem, old_df - 1);
            }
            // else do the generation
            set_x_in_range();
            // create a Spline object for the extended knot sequence
            T_sp bs_obj {
                x_in_range_, surrogate_internal_knots_,
                degree_, surrogate_boundary_knots_
            };
            rmat b_mat { bs_obj.derivative(derivs, true) };
            // remove first and last #degree basis functions
            b_mat = b_mat.cols(degree_, b_mat.n_cols - order_);
            // post-processing
            b_mat = clamp_basis(b_mat);
            if (complete_basis) {
                return b_mat;
            }
            // else
            return mat_wo_col1(b_mat);
        }

        inline rmat integral(
            const bool complete_basis = true
            ) override
        {
            update_knot_sequence();
            set_x_in_range();
            // create a Spline object for the extended knot sequence
            T_sp bs_obj {
                x_in_range_, surrogate_internal_knots_,
                degree_, surrogate_boundary_knots_
            };
            rmat b_mat { bs_obj.integral(true) };
            // remove first and last #degree basis functions
            b_mat = b_mat.cols(degree_, b_mat.n_cols - order_);
            // get initial values at the boundary knots
            rmat v0 {
                bs_obj.set_x(boundary_knots_)->integral(true)
            };
            // remove first and last #degree basis functions
            v0 = v0.cols(degree_, v0.n_cols - order_);
            // clear initial values
            for (size_t i {0}; i < v0.n_cols; ++i) {
                b_mat.col(i) -= v0(0, i);
            }
            // post-processing
            b_mat = clamp_basis(b_mat);
            // get cumulative sum of integral from left boundary knot
            for (size_t j {0}; j < b_mat.n_cols; ++j) {
                b_mat.col(j) = (x_num_shift_ >= 0) %
                    (b_mat.col(j) + v0(1, j) * x_num_shift_);
            }
            // return
            if (complete_basis) {
                return b_mat;
            }
            // else
            return mat_wo_col1(b_mat);
        }

    };                          // end of PeriodicSpline

    // from the template
    using PeriodicBSpline = PeriodicSpline<BSpline>;
    using PeriodicMSpline = PeriodicSpline<MSpline>;

}  // splines2


#endif /* SPLINES2_PERIODICSPLINE_H */
