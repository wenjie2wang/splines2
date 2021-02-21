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

#ifndef SPLINES2_PERIODICMSPLINE_H
#define SPLINES2_PERIODICMSPLINE_H

#include <stdexcept>

#include "common.h"
#include "utils.h"
#include "SplineBase.h"
#include "MSpline.h"

namespace splines2 {

    // define a class for nonnegative periodic spline basis over [a, b]
    // with unit integral over [a, b] based on M-splines
    class PeriodicMSpline
    {
    protected:
        unsigned int degree_ = 3;
        unsigned int order_ = 4;

        // number of internal knots must be >= degree
        rvec internal_knots_;    // (distinct) internal knots
        rvec boundary_knots_;    // [a, b] defining period
        double range_size_;      // b - a
        // [min(knot_sequence), max(knot_sequence)]
        rvec extended_boundary_knots_;
        rvec extended_internal_knots_;

        // internally extended knot sequence
        rvec knot_sequence_;
        bool is_knot_sequence_latest_ = false;

        // number of complete basis: number of internal knots + 1
        unsigned int spline_df_;

        rvec x_;                // input x
        rvec x_in_range_;       // x in range
        rvec x_num_shift_;      // x_ = x_in_range + x_num_shift * range_size
        bool is_x_in_range_latest_ = false;

        // complete spline basis matrix
        rmat spline_basis_;
        bool is_basis_latest_ = false;

        // pre-process some inputs
        // check knots, and do assignment
        inline void clean_knots(const rvec& internal_knots = rvec(),
                                const rvec& boundary_knots = rvec())
        {
            // check before assigment
            if (boundary_knots.has_nan()) {
                throw std::range_error("Boundary knots cannot contain 'NA'.");
            }
            // for unspecified boundary knots
            // 1. do generation if no boundary knots have been set
            // 2. or skip checks if a set of boundary knots have been set
            if (boundary_knots.n_elem == 0) {
                if (boundary_knots_.n_elem != 2 && x_.n_elem > 0) {
                    // set boundary knots to be min(x) and max(x)
                    double left { arma::min(x_) };
                    double right { arma::max(x_) };
                    // check if boundary knots are different
                    if (left == right) {
                        throw std::range_error(
                            "Cannot set boundary knots from 'x'."
                            );
                    }
                    boundary_knots_ = arma::zeros(2);
                    boundary_knots_(0) = left;
                    boundary_knots_(1) = right;
                    range_size_ = boundary_knots_(1) - boundary_knots_(0);
                }
            } else {
                // specified boundary knots
                rvec uni_boundary_knots { arma::unique(boundary_knots) };
                if (uni_boundary_knots.n_elem != 2) {
                    throw std::length_error(
                        "Need two distinct boundary knots.");
                }
                boundary_knots_ = uni_boundary_knots;
                range_size_ = boundary_knots_(1) - boundary_knots_(0);
            }
            if (internal_knots.has_nan()) {
                throw std::range_error(
                    "The internal knots cannot contain 'NA'."
                    );
            }
            // for non-empty internal knots
            // 1. get unique values
            // 2. checks if outside of boundary
            if (internal_knots.n_elem > 0) {
                rvec uni_internal_knots { arma::unique(internal_knots) };
                // check internal knots are inside of boundary knots
                double min_int_knots { uni_internal_knots(0) };
                double max_int_knots {
                    uni_internal_knots(uni_internal_knots.n_elem - 1)
                };
                if (boundary_knots_[0] >= min_int_knots ||
                    boundary_knots_[1] <= max_int_knots) {
                    throw std::range_error(
                        "The internal knots must be inside boundary knots."
                        );
                }
                internal_knots_ = uni_internal_knots;
            } else {
                internal_knots_ = internal_knots;
            }
        }

        // compute spline df
        inline void update_spline_df()
        {
            spline_df_ = internal_knots_.n_elem + 1;
        }

        // extend knot sequence for period spline
        // reference:
        // - Farin, G., & Hansford, D. (2000). The essentials of CAGD.
        // - Schumaker, L. (2007). Spline functions: Basic theory.
        // - splines::splineDesign() when outer.ok = TRUE
        // - Piegl, L., & Tiller, W. (1997). The NURBS book. Springer
        inline void extend_knot_sequence()
        {
            // number of internal knots must >= degree
            if (internal_knots_.n_elem + 1 < degree_) {
                throw std::length_error(
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
            extended_boundary_knots_ = arma::zeros(2);
            extended_boundary_knots_(0) = res(0);
            extended_boundary_knots_(1) = res(res.n_elem - 1);
            extended_internal_knots_ = res.subvec(1, res.n_elem - 2);
            // add multiplicities to "boundary" knots
            knot_sequence_ = arma::zeros(res.n_elem + 2 * degree_);
            for (size_t i {0}; i < knot_sequence_.n_elem; ++i) {
                if (i < order_) {
                    knot_sequence_(i) = extended_boundary_knots_(0);
                } else if (i < knot_sequence_.n_elem - order_) {
                    knot_sequence_(i) = extended_internal_knots_(i - order_);
                } else {
                    knot_sequence_(i) = extended_boundary_knots_(1);
                }
            }
        }

        inline void update_knot_sequence()
        {
            if (! is_knot_sequence_latest_ || knot_sequence_.n_elem == 0) {
                extend_knot_sequence();
                is_knot_sequence_latest_ = true;
            }
        }

        inline void set_x_in_range()
        {
            if (is_x_in_range_latest_) {
                return;
            }
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
        PeriodicMSpline() {}
        virtual ~PeriodicMSpline() {}

        explicit PeriodicMSpline(const SplineBase* pSplineBase)
        {
            x_ = pSplineBase->get_x();
            internal_knots_ = pSplineBase->get_internal_knots();
            boundary_knots_ = pSplineBase->get_boundary_knots();
            degree_ = pSplineBase->get_degree();
            order_ = degree_ + 1;
            update_spline_df();
        }

        // given boundary_knots for consistency with SplineBase
        PeriodicMSpline(const rvec& x,
                        const rvec& internal_knots,
                        const unsigned int degree = 3,
                        const rvec& boundary_knots = rvec())
        {
            x_ = x;
            degree_ = degree;
            clean_knots(internal_knots, boundary_knots);
            if (internal_knots_.n_elem + 1 < degree_) {
                throw std::length_error(
                    "The number of unique internal knots must be >= degree - 1."
                    );
            }
            order_ = degree_ + 1;
            update_spline_df();
        }

        PeriodicMSpline(const rvec& x,
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
            clean_knots(rvec(), boundary_knots);
            // get quantiles of x in range
            set_x_in_range();
            rvec internal_knots { arma_quantile(x_in_range_, prob_vec) };
            clean_knots(internal_knots);
        }


        // function members
        // "setter" functions
        inline PeriodicMSpline* set_x(const rvec& x)
        {
            x_ = x;
            is_x_in_range_latest_ = false;
            is_basis_latest_ = false;
            return this;
        }
        inline PeriodicMSpline* set_x(const double x)
        {
            x_ = num2vec(x);
            is_x_in_range_latest_ = false;
            is_basis_latest_ = false;
            return this;
        }
        inline PeriodicMSpline* set_internal_knots(
            const rvec& internal_knots
            )
        {
            clean_knots(internal_knots);
            // update spline df
            update_spline_df();
            is_knot_sequence_latest_ = false;
            is_basis_latest_ = false;
            return this;
        }
        inline PeriodicMSpline* set_boundary_knots(
            const rvec& boundary_knots
            )
        {
            clean_knots(internal_knots_, boundary_knots);
            is_knot_sequence_latest_ = false;
            is_x_in_range_latest_ = false;
            is_basis_latest_ = false;
            return this;
        }
        inline PeriodicMSpline* set_degree(
            const unsigned int degree
            )
        {
            degree_ = degree;
            order_ = degree + 1;
            is_knot_sequence_latest_ = false;
            is_basis_latest_ = false;
            return this;
        }
        inline PeriodicMSpline* set_order(const unsigned int order)
        {
            if (order > 0) {
                set_degree(order - 1);
            } else {
                throw std::range_error("The 'order' must be at least 1.");
            }
            return this;
        }

        // "getter" functions
        inline rvec get_x() const
        {
            return x_;
        }
        inline rvec get_internal_knots() const
        {
            return internal_knots_;
        }
        inline rvec get_boundary_knots() const
        {
            return boundary_knots_;
        }
        inline unsigned int get_degree() const
        {
            return degree_;
        }
        inline unsigned int get_order() const
        {
            return order_;
        }
        inline unsigned int get_spline_df()
        {
            update_spline_df();
            return spline_df_;
        }


        //! Compute periodic spline basis based on M-spline basis
        //!
        //! @param complete_basis A `bool` value indicating whether to return a
        //! complete spline basis
        //!
        //! @return arma::mat
        inline rmat basis(const bool complete_basis = true)
        {
            // early exit if latest
            if (is_basis_latest_) {
                if (complete_basis) {
                    return spline_basis_;
                }
                // else
                return mat_wo_col1(spline_basis_);
            }
            // else do the generation
            update_knot_sequence();
            set_x_in_range();
            // create a MSpline object for the extended knot sequence
            MSpline ms_obj { x_in_range_, degree_, knot_sequence_ };
            rmat b_mat { ms_obj.basis(true) };
            // remove first and last #degree basis functions
            b_mat = b_mat.cols(degree_, b_mat.n_cols - order_);
            // post-processing
            spline_basis_ = clamp_basis(b_mat);
            is_basis_latest_ = true;
            if (complete_basis) {
                return spline_basis_;
            }
            // else
            return mat_wo_col1(spline_basis_);
        }

        inline rmat derivative(
            const unsigned int derivs = 1,
            const bool complete_basis = true
            )
        {
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
            update_knot_sequence();
            set_x_in_range();
            // create a MSpline object for the extended knot sequence
            MSpline ms_obj {
                x_in_range_, extended_internal_knots_,
                degree_, extended_boundary_knots_
            };
            rmat b_mat { ms_obj.derivative(derivs, true) };
            // remove first and last #degree basis functions
            b_mat = b_mat.cols(degree_, b_mat.n_cols - order_);
            // post-processing
            b_mat = clamp_basis(b_mat);
            is_basis_latest_ = true;
            if (complete_basis) {
                return b_mat;
            }
            // else
            return mat_wo_col1(b_mat);
        }

        inline rmat integral(
            const bool complete_basis = true
            )
        {
            update_knot_sequence();
            set_x_in_range();
            // create a MSpline object for the extended knot sequence
            MSpline ms_obj {
                x_in_range_, extended_internal_knots_,
                degree_, extended_boundary_knots_
            };
            rmat b_mat { ms_obj.integral(true) };
            // remove first and last #degree basis functions
            b_mat = b_mat.cols(degree_, b_mat.n_cols - order_);
            // get initial values at the left boundary knot
            rmat v0 {
                ms_obj.set_x(boundary_knots_(0))->integral(true)
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
                b_mat.col(j) = (b_mat.col(j) + x_num_shift_) %
                    (x_num_shift_ >= 0);
            }
            // return
            if (complete_basis) {
                return b_mat;
            }
            // else
            return mat_wo_col1(b_mat);
        }

    };

}  // splines2


#endif /* SPLINES2_PERIODICMSPLINE_H */
