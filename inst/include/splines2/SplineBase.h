//
// R package splines2 by Wenjie Wang and Jun Yan
// Copyright (C) 2016-2020
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

#ifndef SPLINE_BASE_H
#define SPLINE_BASE_H

#include <stdexcept>

#include "common.h"

namespace splines2 {

    // define base class for some regression splines
    class SplineBase {
    protected:
        rvec x;
        ivec x_index;
        rvec internal_knots;
        rvec boundary_knots;
        rvec knot_sequence;
        unsigned int degree = 3;
        unsigned int order = 4;
        unsigned int df = 0;    // degree of freedom

        // Does basis need updating
        bool is_basis_latest = false;
        // Do the intermediate results need updating
        bool is_knot_sequence_latest = false;
        bool is_x_index_latest = false;

        rmat full_basis;        // complete spline matrix

        // pre-process some inputs
        inline void clean_internal_knots()
        {
            this->internal_knots = Rcpp::sort_unique(this->internal_knots);
        }
        inline void clean_boundary_knots()
        {
            this->boundary_knots = Rcpp::sort_unique(this->boundary_knots);
            if (this->boundary_knots.size() > 2) {
                this->boundary_knots = Rcpp::head(this->boundary_knots, 2);
                // FIXME: throw warning?
            }
        }
        inline rvec get_knot_sequence(const unsigned int order = 1)
        {
            rvec out { rvec(internal_knots.size() + 2 * order) };
            for (long i {0}; i < out.size(); ++i) {
                if (i < order) {
                    out(i) = boundary_knots(0);
                } else if (i < out.size() - order) {
                    out(i) = internal_knots(i - order);
                } else {
                    out(i) = boundary_knots(1);
                }
            }
            return out;
        }
        inline void update_knot_sequence()
        {
            if (! is_knot_sequence_latest || knot_sequence.size() == 0) {
                this->knot_sequence = get_knot_sequence(this->order);
            }
            is_knot_sequence_latest = true;
        }
        inline void update_x_index()
        {
            if (! is_x_index_latest || x_index.size() == 0) {
                x_index = ivec(x.size());
                rvec uni_knots { get_knot_sequence() };
                for (long i {0}; i < x.size(); ++i) {
                    // assume x all inside of range(knots)
                    long left_index {0};
                    long right_index { uni_knots.size() - 1 };
                    while (right_index > left_index + 1) {
                        long cur_index {
                            (left_index + right_index) / 2
                        };
                        if (x(i) >= uni_knots(cur_index)) {
                            left_index = cur_index;
                        } else {
                            right_index = cur_index;
                        }
                    }
                    // close on right boundary knot
                    x_index(i) = left_index;
                }
            }
            is_x_index_latest = true;
        }

    public:
        // the default constructor
        SplineBase() {}

        // main constructor
        SplineBase(const rvec& x,
                   const unsigned int degree = 3,
                   const rvec& internal_knots = rvec(),
                   const rvec& boundary_knots = rvec()) :
            x { x },
            internal_knots { internal_knots },
            boundary_knots { boundary_knots },
            degree { degree }
        {
            // by default, no internal knots
            if (this->boundary_knots.size() == 0) {
                // set boundary knots to be min(x) and max(x)
                this->boundary_knots = rvec(2);
                this->boundary_knots(0) = Rcpp::min(x);
                this->boundary_knots(1) = Rcpp::max(x);
            } else {
                clean_boundary_knots();
            }
            this->order = this->degree + 1;
            if (internal_knots.size() == 0) {
                this->df = this->order;
            } else {
                this->df = internal_knots.size() + this->order;
                clean_internal_knots();
            }
        }

        // explicit constructor
        explicit SplineBase(const rvec& x)
        {
            SplineBase(x, 3, rvec(), rvec());
        }

        // function members
        // the "setter" functions
        inline void set_x(const rvec& x)
        {
            this->x = x;
            this->is_x_index_latest = false;
            this->is_basis_latest = false;
        }
        inline void set_internal_knots(const rvec& internal_knots)
        {
            this->internal_knots = internal_knots;
            clean_internal_knots();
            this->is_knot_sequence_latest = false;
            this->is_x_index_latest = false;
            this->is_basis_latest = false;
        }
        inline void set_boundary_knots(const rvec& boundary_knots)
        {
            this->boundary_knots = boundary_knots;
            clean_boundary_knots();
            this->is_knot_sequence_latest = false;
            this->is_x_index_latest = false;
            this->is_basis_latest = false;
        }
        inline void set_degree(const unsigned int degree)
        {
            this->degree = degree;
            this->order = degree + 1;
            this->is_knot_sequence_latest = false;
            this->is_basis_latest = false;
        }
        inline void set_order(const unsigned int order)
        {
            if (order > 0) {
                this->set_degree(order - 1);
            } else {
                throw std::range_error("The 'order' must be at least 1.");
            }
        }

        // the "getter" functions
        inline rvec get_x() const
        {
            return this->x;
        }
        inline rvec get_internal_knots() const
        {
            return this->internal_knots;
        }
        inline rvec get_boundary_knots() const
        {
            return this->boundary_knots;
        }
        inline unsigned int get_degree() const
        {
            return this->degree;
        }
        inline unsigned int get_order() const
        {
            return this->order;
        }
        inline unsigned int get_df() const
        {
            return this->df;
        }


    };

}

#endif
