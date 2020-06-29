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
#include "utils.h"

namespace splines2arma {

    // define base class for some regression splines
    class SplineBase {
    protected:
        // set and get
        rvec x_;
        rvec internal_knots_;
        rvec boundary_knots_;
        unsigned int degree_ = 3;
        unsigned int order_ = 4;
        // degree of freedom of complete spline basis
        // notice that this argument is different with the df in splines::bs()
        unsigned int spline_df_ = 4;

        // knot sequence
        rvec knot_sequence_;
        bool is_knot_sequence_latest_ = false;

        // index of x relative to internal knots
        uvec x_index_;
        bool is_x_index_latest_ = false;

        // complete spline matrix
        rmat spline_basis_;
        bool is_basis_latest_ = false;


        // pre-process some inputs
        // check knots, and do assignment
        inline void clean_knots(const arma::vec& internal_knots = nan_vec(),
                                const arma::vec& boundary_knots = rvec())
        {
            if (boundary_knots.n_elem == 0 && this->x_.n_elem > 0 &&
                this->boundary_knots_.n_elem != 2) {
                // set boundary knots to be min(x) and max(x)
                this->boundary_knots_ = arma::zeros(2);
                this->boundary_knots_(0) = arma::min(this->x_);
                this->boundary_knots_(1) = arma::max(this->x_);
                // check if boundary knots are different
                if (this->boundary_knots_(0) == this->boundary_knots_(1)) {
                    throw std::range_error(
                        "Cannot set boundary knots from unique x."
                        );
                }
            } else {
                arma::vec uni_boundary_knots { arma::unique(boundary_knots) };
                if (uni_boundary_knots.n_elem != 2) {
                    throw std::length_error(
                        "Need two distinct boundary knots."
                        );
                }
                this->boundary_knots_ = uni_boundary_knots;
            }
            if (internal_knots.has_nan()) {
                // do nothing
            } else if (internal_knots.n_elem > 0) {
                rvec uni_internal_knots { arma::unique(internal_knots) };
                // check internal knots are inside of boundary knots
                double min_int_knots { arma::min(internal_knots) };
                double max_int_knots { arma::max(internal_knots) };
                if (this->boundary_knots_[0] >= min_int_knots ||
                    this->boundary_knots_[1] <= max_int_knots) {
                    throw std::range_error(
                        "Internal knots must be set inside of boundary knots."
                        );
                }
                this->internal_knots_ = uni_internal_knots;
            } else {
                this->internal_knots_ = internal_knots;
            }
        }

        // compute spline df
        inline void update_spline_df()
        {
            this->spline_df_ = this->internal_knots_.n_elem + this->order_;
        }

        inline arma::vec get_knot_sequence(const unsigned int order = 1)
        {
            arma::vec out { arma::zeros(internal_knots_.n_elem + 2 * order) };
            for (size_t i {0}; i < out.n_elem; ++i) {
                if (i < order) {
                    out(i) = boundary_knots_(0);
                } else if (i < out.n_elem - order) {
                    out(i) = internal_knots_(i - order);
                } else {
                    out(i) = boundary_knots_(1);
                }
            }
            return out;
        }
        inline void update_knot_sequence()
        {
            if (! is_knot_sequence_latest_ || knot_sequence_.n_elem == 0) {
                this->knot_sequence_ = get_knot_sequence(this->order_);
                is_knot_sequence_latest_ = true;
            }
        }
        // allow x outside of boundary
        // let degree 0 basis take 1 outside boundary
        inline void update_x_index()
        {
            if (! is_x_index_latest_ || x_index_.n_elem == 0) {
                x_index_ = arma::zeros<arma::uvec>(x_.n_elem);
                for (size_t i {0}; i < x_.n_elem; ++i) {
                    size_t left_index {0};
                    size_t right_index { internal_knots_.n_elem };
                    while (right_index > left_index) {
                        size_t cur_index {
                            (left_index + right_index) / 2
                        };
                        if (x_(i) < internal_knots_(cur_index)) {
                            right_index = cur_index;
                        } else {
                            left_index = cur_index + 1;
                        }
                    }
                    x_index_(i) = left_index;
                }
                is_x_index_latest_ = true;
            }
        }

    public:
        // the default constructor
        SplineBase() {}

        // explicit constructor
        explicit SplineBase(const arma::vec& x) :
            SplineBase(x, arma::vec(), 3, arma::vec()) {}

        // constructor with specificied internal_knots
        SplineBase(const arma::vec& x,
                   const arma::vec& internal_knots = arma::vec(),
                   const unsigned int degree = 3,
                   const arma::vec& boundary_knots = arma::vec()) :
            x_ { x },
            degree_ { degree }
        {
            this->clean_knots(internal_knots, boundary_knots);
            this->order_ = this->degree_ + 1;
            this->update_spline_df();
        }

        // constructor with specified df
        SplineBase(const arma::vec& x,
                   const unsigned int spline_df = 4,
                   const unsigned int degree = 3,
                   const arma::vec& boundary_knots = arma::vec()) :
            x_ { x },
            degree_ { degree }
        {
            this->order_ = this->degree_ + 1;
            if (spline_df < this->order_) {
                throw std::range_error("The specified df was too small.");
            }
            this->spline_df_ = spline_df;
            // determine internal knots by spline_df and x
            unsigned int n_internal_knots { this->spline_df_ - this->order_ };
            if (n_internal_knots == 0) {
                this->clean_knots(rvec(), boundary_knots);
            } else {
                rvec prob_vec { arma::linspace(0, 1, n_internal_knots + 2) };
                prob_vec = prob_vec.subvec(1, n_internal_knots);
                rvec internal_knots { arma_quantile(this->x_, prob_vec) };
                this->clean_knots(rvec(), boundary_knots);
            }
        }

        // function members
        // "setter" functions
        inline SplineBase* set_x(const arma::vec& x)
        {
            this->x_ = x;
            this->is_x_index_latest_ = false;
            this->is_basis_latest_ = false;
            return this;
        }
        inline SplineBase* set_internal_knots(const arma::vec& internal_knots)
        {
            clean_knots(internal_knots);
            // update spline df
            this->update_spline_df();
            this->is_knot_sequence_latest_ = false;
            this->is_x_index_latest_ = false;
            this->is_basis_latest_ = false;
            return this;
        }
        inline SplineBase* set_boundary_knots(const arma::vec& boundary_knots)
        {
            clean_knots(nan_vec(), boundary_knots);
            this->is_knot_sequence_latest_ = false;
            this->is_x_index_latest_ = false;
            this->is_basis_latest_ = false;
            return this;
        }
        inline SplineBase* set_degree(const unsigned int degree)
        {
            this->degree_ = degree;
            this->order_ = degree + 1;
            // update spline df
            this->update_spline_df();
            this->is_knot_sequence_latest_ = false;
            this->is_basis_latest_ = false;
            return this;
        }
        inline SplineBase* set_order(const unsigned int order)
        {
            if (order > 0) {
                this->set_degree(order - 1);
            } else {
                throw std::range_error("The 'order' must be at least 1.");
            }
            return this;
        }

        // "getter" functions
        inline arma::vec get_x() const
        {
            return this->x_;
        }
        inline arma::vec get_internal_knots() const
        {
            return this->internal_knots_;
        }
        inline arma::vec get_boundary_knots() const
        {
            return this->boundary_knots_;
        }
        inline unsigned int get_degree() const
        {
            return this->degree_;
        }
        inline unsigned int get_order() const
        {
            return this->order_;
        }
        inline unsigned int get_spline_df() const
        {
            return this->spline_df_;
        }


    };

}

#endif
