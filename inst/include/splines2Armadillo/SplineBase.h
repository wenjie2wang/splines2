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

#ifndef SPLINE2_SPLINE_BASE_H
#define SPLINE2_SPLINE_BASE_H

#include <stdexcept>
#include <algorithm>

#include "common.h"
#include "utils.h"

namespace splines2 {

    // define base class for some regression splines
    class SplineBase
    {
    protected:
        // setter and getter
        rvec x_ = rvec();
        rvec internal_knots_ = rvec();
        rvec boundary_knots_ = rvec();
        unsigned int degree_ = 3;
        unsigned int order_ = 4;

        // degree of freedom of complete spline basis
        // notice that this argument is different with the df in splines::bs()
        unsigned int spline_df_ = 4;

        // knot sequence
        rvec knot_sequence_ = rvec();
        bool has_internal_multiplicity_ = false;
        bool is_knot_sequence_latest_ = false;
        bool is_extended_knot_sequence_ = false;

        // for extended knot sequence
        // [min(knot_sequence), max(knot_sequence)]
        rvec surrogate_internal_knots_;
        rvec surrogate_boundary_knots_;

        // index of x relative to internal knots
        uvec x_index_ = uvec();
        bool is_x_index_latest_ = false;

        // pre-process some inputs
        // check knots, and do assignment
        inline virtual void simplify_knots(
            const rvec& internal_knots = rvec(),
            const rvec& boundary_knots = rvec()
            )
        {
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
                            "Cannot set boundary knots from x."
                            );
                    }
                    boundary_knots_ = arma::zeros(2);
                    boundary_knots_(0) = left;
                    boundary_knots_(1) = right;
                }
            } else {
                // check before assigment
                if (boundary_knots.has_nan()) {
                    throw std::range_error("Boundary knots cannot contain NA.");
                }
                // specified boundary knots
                rvec uni_boundary_knots { arma::unique(boundary_knots) };
                if (uni_boundary_knots.n_elem != 2) {
                    throw std::range_error(
                        "Need two distinct boundary knots.");
                }
                boundary_knots_ = uni_boundary_knots;
            }
            if (internal_knots.has_nan()) {
                throw std::range_error("Internal knots cannot contain NA.");
            }
            // for non-empty internal knots
            if (internal_knots.n_elem > 0) {
                // check internal knots are inside boundary knots
                rvec sorted_internal_knots { arma::sort(internal_knots) };
                double min_int_knots { sorted_internal_knots(0) };
                double max_int_knots {
                    sorted_internal_knots(sorted_internal_knots.n_elem - 1)
                };
                if (boundary_knots_.n_elem == 2 &&
                    (boundary_knots_[0] >= min_int_knots ||
                     boundary_knots_[1] <= max_int_knots)) {
                    throw std::range_error(
                        "Internal knots must be set inside boundary."
                        );
                }
                // check multiplicity
                rvec tmp {
                    arma::join_vert(sorted_internal_knots, boundary_knots_)
                };
                has_internal_multiplicity_ = any_duplicated(tmp);
                internal_knots_ = sorted_internal_knots;
            } else {
                has_internal_multiplicity_ = false;
                internal_knots_ = rvec();
            }
        }

        // compute spline df
        inline virtual void update_spline_df()
        {
            spline_df_ = internal_knots_.n_elem + order_;
        }

        // get simple knot sequence
        inline virtual rvec get_simple_knot_sequence(
            const rvec& internal_knots,
            const rvec& boundary_knots,
            const unsigned int order
            ) const
        {
            rvec out { arma::zeros(internal_knots.n_elem + 2 * order) };
            rvec::iterator it { out.begin() }, it_end { out.end() - 1 };
            rvec::const_iterator ii { internal_knots.begin() };
            for (size_t i {0}; i < order; ++i, ++it, --it_end) {
                *it = boundary_knots(0);
                *it_end = boundary_knots(1);
            }
            for (++it_end; it != it_end; ++it, ++ii) {
                *it = *ii;
            }
            return out;
        }

        // set simple knot sequence
        inline virtual void set_simple_knot_sequence()
        {
            knot_sequence_ = get_simple_knot_sequence(
                internal_knots_, boundary_knots_, order_
                );
            is_knot_sequence_latest_ = true;
        }

        // set extended knot sequence
        inline virtual void set_extended_knot_sequence(
            const rvec& knot_sequence
            )
        {
            // check the length of specified knot sequence
            if (knot_sequence.n_elem < 2 * order_) {
                throw std::range_error(
                   "The length of specified knot sequence is too small."
                    );
            }
            unsigned int n_internal_knots {
                knot_sequence.n_elem - 2 * order_
            };
            // sort knot sequence
            knot_sequence_ = arma::sort(knot_sequence);
            // get boundary knots
            boundary_knots_ = arma::zeros(2);
            boundary_knots_(0) = knot_sequence_(degree_);
            boundary_knots_(1) = knot_sequence_(knot_sequence_.n_elem - order_);
            if (isAlmostEqual(boundary_knots_(0), boundary_knots_(1))) {
                throw std::range_error(
                    "The specified knot sequence has the same boundary knots."
                    );
            }
            // get internal knots
            if (n_internal_knots > 0) {
                internal_knots_ = knot_sequence_.subvec(
                    order_, order_ + n_internal_knots - 1);
                // check multiplicity
                rvec tmp {
                    arma::join_vert(internal_knots_, boundary_knots_)
                };
                has_internal_multiplicity_ = any_duplicated(tmp);
            } else {
                internal_knots_ = rvec();
                has_internal_multiplicity_ = false;
            }
            // set surrogate knots
            surrogate_boundary_knots_ = arma::zeros(2);
            surrogate_boundary_knots_(0) = knot_sequence_(0);
            surrogate_boundary_knots_(1) =
                knot_sequence_(knot_sequence_.n_elem - 1);
            surrogate_internal_knots_ =
                knot_sequence_.subvec(1, knot_sequence_.n_elem - 2);
            // check if it is actually a simple knot sequence
            is_extended_knot_sequence_ = ! (
                isAlmostEqual(boundary_knots_(0),
                              surrogate_boundary_knots_(0)) &&
                isAlmostEqual(boundary_knots_(1),
                              surrogate_boundary_knots_(1))
                ) || has_internal_multiplicity_;
            // set flags to prevent knot sequence from being updated
            is_knot_sequence_latest_ = true;
        }

        inline virtual void update_knot_sequence()
        {
            if (is_knot_sequence_latest_ && knot_sequence_.n_elem > 0) {
                return;
            }
            if (is_extended_knot_sequence_) {
                set_extended_knot_sequence(knot_sequence_);
            } else {
                set_simple_knot_sequence();
            }
        }

        // allow x outside of boundary
        // let degree 0 basis take 1 outside boundary
        // the way to generate knot index for each x
        // avoids denominator = 0
        inline virtual void update_x_index()
        {
            if (is_x_index_latest_ && x_index_.n_elem > 0) {
                return;
            }
            x_index_ = arma::zeros<arma::uvec>(x_.n_elem);
            arma::uvec::iterator xi { x_index_.begin() };
            arma::vec::iterator it { x_.begin() }, pos, ie { x_.end() },
                knots_begin { internal_knots_.begin() },
                knots_end { internal_knots_.end() };
            for (; it != ie; ++it, ++xi) {
                pos = std::upper_bound(knots_begin, knots_end, *it);
                *xi = std::distance(knots_begin, pos);
            }
            is_x_index_latest_ = true;
        }

        // check if simple knot sequence
        inline virtual void stopifnot_simple_knot_sequence() const
        {
            if (has_internal_multiplicity_ || is_extended_knot_sequence_) {
                throw std::range_error(
                    "Expected a simple knot sequence."
                    );
            }
        }

    public:
        // the default constructor
        SplineBase() {}
        virtual ~SplineBase() {}

        // explicit constructor
        explicit SplineBase(const SplineBase* pSplineBase) :
            x_ { pSplineBase->x_ },
            internal_knots_ { pSplineBase->internal_knots_ },
            boundary_knots_ { pSplineBase->boundary_knots_ },
            degree_ { pSplineBase->degree_  },
            knot_sequence_ { pSplineBase->knot_sequence_ },
            has_internal_multiplicity_ {
                pSplineBase->has_internal_multiplicity_ },
            is_knot_sequence_latest_ {
                pSplineBase->is_knot_sequence_latest_ },
            is_extended_knot_sequence_ {
                pSplineBase->is_extended_knot_sequence_ },
            surrogate_internal_knots_ {
                pSplineBase->surrogate_internal_knots_ },
            surrogate_boundary_knots_ {
                pSplineBase->surrogate_boundary_knots_ },
            x_index_ { pSplineBase->x_index_ },
            is_x_index_latest_ { pSplineBase->is_x_index_latest_ }
        {
            order_ = degree_ + 1;
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

        // constructor with specificied internal_knots
        SplineBase(const rvec& x,
                   const rvec& internal_knots,
                   const unsigned int degree = 3,
                   const rvec& boundary_knots = rvec()) :
            x_ { x },
            degree_ { degree }
        {
            simplify_knots(internal_knots, boundary_knots);
            order_ = degree_ + 1;
        }

        // constructor with specified df
        SplineBase(const rvec& x,
                   const unsigned int spline_df,
                   const unsigned int degree = 3,
                   const rvec& boundary_knots = rvec()) :
            x_ { x },
            degree_ { degree }
        {
            order_ = degree_ + 1;
            if (spline_df < order_) {
                throw std::range_error("The specified df was too small.");
            }
            spline_df_ = spline_df;
            // determine internal knots by spline_df and x
            unsigned int n_internal_knots { spline_df_ - order_ };
            if (n_internal_knots == 0) {
                simplify_knots(rvec(), boundary_knots);
            } else {
                rvec prob_vec { arma::linspace(0, 1, n_internal_knots + 2) };
                prob_vec = prob_vec.subvec(1, n_internal_knots);
                simplify_knots(rvec(), boundary_knots);
                // get quantiles of x within boundary only
                rvec x_inside { get_inside_x(x, boundary_knots_) };
                rvec internal_knots { quantile(x_inside, prob_vec) };
                simplify_knots(internal_knots);
            }
        }

        // constructor for a given knot sequence
        SplineBase(const rvec& x,
                   const unsigned int degree,
                   const rvec& knot_sequence)
        {
            x_ = x;
            degree_ = degree;
            order_ = degree_ + 1;
            set_extended_knot_sequence(knot_sequence);
        }

        // function members
        // "setter" functions
        inline virtual SplineBase* set_x(const rvec& x)
        {
            x_ = x;
            is_x_index_latest_ = false;
            return this;
        }
        inline virtual SplineBase* set_x(const double x)
        {
            x_ = num2vec(x);
            is_x_index_latest_ = false;
            return this;
        }
        inline virtual SplineBase* set_internal_knots(
            const rvec& internal_knots
            )
        {
            if (! is_approx_equal(internal_knots_, internal_knots)) {
                simplify_knots(internal_knots);
                update_spline_df();
                is_knot_sequence_latest_ = false;
                is_x_index_latest_ = false;
            }
            return this;
        }
        inline virtual SplineBase* set_boundary_knots(
            const rvec& boundary_knots
            )
        {
            if (! is_approx_equal(boundary_knots_, boundary_knots)) {
                simplify_knots(internal_knots_, boundary_knots);
                is_knot_sequence_latest_ = false;
                is_x_index_latest_ = false;
            }
            return this;
        }
        inline virtual SplineBase* set_knot_sequence(
            const rvec& knot_sequence
            )
        {
            if (! is_approx_equal(knot_sequence_, knot_sequence)) {
                set_extended_knot_sequence(knot_sequence);
            }
            return this;
        }
        inline virtual SplineBase* set_degree(
            const unsigned int degree
            )
        {
            if (degree_ != degree) {
                degree_ = degree;
                order_ = degree + 1;
                update_spline_df();
                if (is_extended_knot_sequence_) {
                    // if a knot sequence has been set and it is extended,
                    // update knot sequence
                    // so that internal/boundary knots are update-to-date
                    set_extended_knot_sequence(knot_sequence_);
                } else {
                    is_knot_sequence_latest_ = false;
                }
            }
            return this;
        }
        inline virtual SplineBase* set_order(const unsigned int order)
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
        inline rvec get_knot_sequence()
        {
            update_knot_sequence();
            return knot_sequence_;
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

        // define pure virtual functions
        inline virtual rmat basis(
            const bool complete_basis = true
            ) = 0;
        inline virtual rmat derivative(
            const unsigned int derivs = 1,
            const bool complete_basis = true
            ) = 0;
        inline virtual rmat integral(
            const bool complete_basis = true
            ) = 0;

    };

}

#endif
