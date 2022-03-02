//
// R package splines2 by Wenjie Wang and Jun Yan
// Copyright (C) 2016-2022
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

#ifndef SPLINES2_BERNSTEINPOLY_H
#define SPLINES2_BERNSTEINPOLY_H

#include <stdexcept>

#include "common.h"
#include "SplineBase.h"
#include "utils.h"

namespace splines2 {

    // define class for generalized Bernstein polynomials over [a, b]
    class BernsteinPoly
    {
    protected:
        unsigned int degree_ = 3;
        unsigned int order_ = 4;
        rvec boundary_knots_;   // [a, b]
        double range_size_ = 1; // b - a
        rvec x_;

        // check x
        inline void check_x(const rvec& x) {
            if (x.has_nan()) {
                throw std::range_error("x cannot contain NA.");
            }
            if (boundary_knots_.n_elem == 2) {
                for (size_t i {0}; i < x.n_elem; ++i) {
                    if (x(i) < boundary_knots_(0) ||
                        x(i) > boundary_knots_(1)) {
                        throw std::range_error(
                            "The 'x' must be inside of boundary."
                            );
                    }
                }
            }
            x_ = x;
        }
        // check boundary
        inline void check_boundary(const rvec& boundary_knots)
        {
            if (boundary_knots.has_nan()) {
                throw std::range_error("Boundary knots cannot contain NA.");
            }
            double left { boundary_knots(0) };
            double right { boundary_knots(1) };
            if (left >= right) {
                throw std::range_error(
                    "The left boundary must be less than the right boundary."
                    );
            }
            boundary_knots_ = arma::zeros(2);
            boundary_knots_(0) = left;
            boundary_knots_(1) = right;
            range_size_ = right - left;
        }
        inline void autoset_x_and_boundary(const rvec& x)
        {
            if (x.n_elem == 0) {
                return;
            } else if (x.has_nan()) {
                throw std::range_error("x cannot contain NA.");
            }
            boundary_knots_ = arma::zeros(2);
            boundary_knots_(0) = arma::min(x);
            boundary_knots_(1) = arma::max(x);
            range_size_ = boundary_knots_(1) - boundary_knots_(0);
            x_ = x;
        }

    public:
        BernsteinPoly() {}
        virtual ~BernsteinPoly() {}

        explicit BernsteinPoly(const BernsteinPoly* pBernsteinPoly)
        {
            x_ = pBernsteinPoly->x_;
            degree_ = pBernsteinPoly->degree_;
            order_ = degree_ + 1;
            if (pBernsteinPoly->boundary_knots_.n_elem == 0) {
                autoset_x_and_boundary(x_);
            } else {
                check_boundary(pBernsteinPoly->boundary_knots_);
            }
        }

        explicit BernsteinPoly(const SplineBase* pSplineBase)
        {
            x_ = pSplineBase->get_x();
            degree_ = pSplineBase->get_degree();
            order_ = degree_ + 1;
            rvec bound_knots { pSplineBase->get_boundary_knots() };
            if (bound_knots.n_elem == 0) {
                autoset_x_and_boundary(x_);
            } else {
                check_boundary(bound_knots);
            }
        }

        // explicit conversion
        template <typename T>
        explicit operator T() const {
            T obj;
            obj.set_x(x_)->
                set_degree(degree_)->
                set_boundary_knots(boundary_knots_);
            return obj;
        }

        // given boundary_knots for consistency with SplineBase
        BernsteinPoly(const rvec& x,
                      const unsigned int degree,
                      const rvec& boundary_knots = rvec()) :
            degree_ { degree },
            order_ { degree + 1 }
        {
            if (boundary_knots.n_elem == 0) {
                autoset_x_and_boundary(x);
            } else if (boundary_knots.n_elem != 2) {
                throw std::range_error("Need two distinct boundary knots.");
            } else {
                check_boundary(boundary_knots);
                check_x(x);
            }
        }

        // setter functions
        inline BernsteinPoly* set_x(const rvec& x)
        {
            check_x(x);
            return this;
        }
        inline BernsteinPoly* set_x(const double x)
        {
            check_x(num2vec(x));
            return this;
        }
        inline BernsteinPoly* set_degree(const unsigned int degree)
        {
            degree_ = degree;
            order_ = degree + 1;
            return this;
        }
        inline BernsteinPoly* set_order(const unsigned int order)
        {
            if (order > 0) {
                set_degree(order - 1);
            } else {
                throw std::range_error("The 'order' must be at least 1.");
            }
            return this;
        }
        // placeholder for conversion
        inline BernsteinPoly* set_internal_knots(const rvec& internal_knots)
        {
            if (internal_knots.n_elem > 0) {
                // do nothing
            }
            return this;
        }
        inline BernsteinPoly* set_boundary_knots(const rvec& boundary_knots)
        {
            check_boundary(boundary_knots);
            check_x(x_);
            return this;
        }

        // getter functions
        inline rvec get_x() const
        {
            return x_;
        }
        inline unsigned int get_degree() const
        {
            return degree_;
        }
        inline unsigned int get_order() const
        {
            return order_;
        }
        inline rvec get_boundary_knots() const
        {
            return boundary_knots_;
        }

        // construct polynomial basis by recursive formula
        inline rmat basis(const bool complete_basis = true)
        {
            // define output matrix
            rmat b_mat {
                arma::ones(x_.n_elem, order_)
                    };
            // only do if degree >= 1
            for (unsigned int k {1}; k <= degree_; ++k) {
                for (size_t i {0}; i < x_.n_elem; ++i) {
                    double saved { 0 };
                    for (size_t j {0}; j < k; ++j) {
                        double term { b_mat(i, j) / range_size_ };
                        b_mat(i, j) = saved +
                            (boundary_knots_(1) - x_(i)) * term;
                        saved = (x_(i) - boundary_knots_(0)) * term;
                    }
                    b_mat(i, k) = saved;
                }
            }
            if (complete_basis) {
                return b_mat;
            }
            // else
            return mat_wo_col1(b_mat);
        }

        // derivatives
        inline rmat derivative(const unsigned int derivs = 1,
                               const bool complete_basis = true)
        {
            if (derivs == 0) {
                throw std::range_error(
                    "'derivs' has to be a positive integer."
                    );
            }
            // early exit if derivs is large enough
            if (degree_ < derivs) {
                if (complete_basis) {
                    return arma::zeros(x_.n_elem, order_);
                }
                if (order_ == 1) {
                    throw std::range_error("No column left in the matrix.");
                }
                return arma::zeros(x_.n_elem, order_ - 1);
            }
            BernsteinPoly bp_obj2 { this };
            // get basis matrix for (degree - derivs)
            bp_obj2.set_degree(degree_ - derivs);
            rmat d_mat { bp_obj2.basis(true) };
            // add zero columns
            d_mat = add_zero_cols(d_mat, order_ - d_mat.n_cols);
            // derivatives by recursive formula
            for (unsigned int k {1}; k <= derivs; ++k) {
                const unsigned int k_offset { derivs - k };
                const size_t numer { degree_ - k_offset };
                const double numer2 { numer / range_size_  };
                for (size_t i {0}; i < x_.n_elem; ++i) {
                    double saved { 0 };
                    for (size_t j {0}; j < numer; ++j) {
                        double term { numer2 * d_mat(i, j) };
                        d_mat(i, j) = saved - term;
                        saved = term;
                    }
                    d_mat(i, numer) = saved;
                }
            }
            // remove the first column if needed
            if (complete_basis) {
                return d_mat;
            }
            // else
            return mat_wo_col1(d_mat);
        }

        // integrals
        inline rmat integral(const bool complete_basis = true)
        {
            BernsteinPoly bp_obj2 { this };
            // get basis matrix for (degree + 1) with intercept
            bp_obj2.set_degree(order_);
            rmat i_mat { bp_obj2.basis(false) };
            // integral by recursive formula
            const double fac { range_size_ / order_ };
            for (unsigned int i {0}; i < x_.n_elem; ++i) {
                arma::rowvec tmp { i_mat.row(i) * fac };
                i_mat.row(i) = rev_cum_sum(tmp) ;
            }
            // remove the first column if needed
            if (complete_basis) {
                return i_mat;
            }
            return mat_wo_col1(i_mat);
        }

    };

}  // splines2


#endif /* SPLINES2_BERNSTEINPOLY_H */
