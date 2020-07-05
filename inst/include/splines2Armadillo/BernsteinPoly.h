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

#ifndef SPLINES2_BERNSTEINPOLY_H
#define SPLINES2_BERNSTEINPOLY_H

#include <stdexcept>

#include "common.h"
#include "utils.h"

namespace splines2 {

    // define class for generalized Bernstein polynomials over [a, b]
    class BernsteinPoly
    {
    protected:
        unsigned int degree_ = 3;
        unsigned int order_ = 4;
        double left_boundary_ = 0;  // a
        double right_boundary_ = 1; // b
        double range_size_ = 1;     // b - a
        rvec x_;

        rmat poly_basis_;
        bool is_basis_latest_ = false;

        // check x
        inline void check_x(const rvec& x) {
            for (size_t i {0}; i < x.n_elem; ++i) {
                if (x(i) < left_boundary_ || x(i) > right_boundary_) {
                    throw std::range_error(
                        "The 'x' must be inside of boundary."
                        );
                }
            }
            this->x_ = x;
        }
        // check boundary
        inline void check_boundary(const double left, const double right)
        {
            if (left >= right) {
                throw std::range_error(
                    "The left boundary must be less than the right boundary."
                    );
            }
            this->left_boundary_ = left;
            this->right_boundary_ = right;
            this->range_size_ = right - left;
        }

    public:
        BernsteinPoly() {}
        virtual ~BernsteinPoly() {}

        explicit BernsteinPoly(const rvec& x)
        {
            this->check_x(x);
        }

        BernsteinPoly(const rvec& x,
                      const unsigned int degree,
                      const double left_boundary = 0,
                      const double right_boundary = 1) :
            degree_ { degree },
            order_ { degree + 1 }
        {
            this->check_boundary(left_boundary, right_boundary);
            this->check_x(x);
        }

        // setter functions
        inline BernsteinPoly* set_x(const rvec& x)
        {
            this->check_x(x);
            this->is_basis_latest_ = false;
            return this;
        }
        inline BernsteinPoly* set_x(const double x)
        {
            this->check_x(num2vec(x));
            this->is_basis_latest_ = false;
            return this;
        }
        inline BernsteinPoly* set_degree(const unsigned int degree)
        {
            this->degree_ = degree;
            this->order_ = degree + 1;
            this->is_basis_latest_ = false;
            return this;
        }
        inline BernsteinPoly* set_order(const unsigned int order)
        {
            if (order > 0) {
                this->set_degree(order - 1);
            } else {
                throw std::range_error("The 'order' must be at least 1.");
            }
            return this;
        }
        inline BernsteinPoly* set_boundary(const double left = 0,
                                           const double right = 1)
        {
            this->check_boundary(left, right);
            return this;
        }

        // getter functions
        inline rvec get_x() const
        {
            return this->x_;
        }
        inline unsigned int get_degree() const
        {
            return this->degree_;
        }
        inline unsigned int get_order() const
        {
            return this->order_;
        }
        inline unsigned int get_left_boundary() const
        {
            return this->left_boundary_;
        }
        inline unsigned int get_right_boundary() const
        {
            return this->right_boundary_;
        }

        // construct polynomial bases by recursive formula
        inline virtual rmat basis(const bool complete_basis = true)
        {
            // early exit if latest
            if (this->is_basis_latest_) {
                if (complete_basis) {
                    return this->poly_basis_;
                }
                // else
                return mat_wo_col1(this->poly_basis_);
            }
            // define output matrix
            rmat b_mat {
                arma::ones(this->x_.n_elem, this->order_)
            };
            // only do if degree >= 1
            for (unsigned int k {1}; k <= degree_; ++k) {
                for (size_t i {0}; i < x_.n_elem; ++i) {
                    double saved { 0 };
                    for (size_t j {0}; j < k; ++j) {
                        double term { b_mat(i, j) / range_size_ };
                        b_mat(i, j) = saved + (right_boundary_ - x_(i)) * term;
                        saved = (x_(i) - left_boundary_) * term;
                    }
                    b_mat(i, k) = saved;
                }
            }
            // prepare to return
            this->poly_basis_ = b_mat;
            this->is_basis_latest_ = true;
            if (complete_basis) {
                return b_mat;
            }
            // else
            return mat_wo_col1(b_mat);
        }

        // derivatives
        inline virtual rmat derivative(const unsigned int derivs = 1,
                                       const bool complete_basis = true)
        {
            if (derivs == 0) {
                throw std::range_error(
                    "'derivs' has to be a positive integer."
                    );
            }
            // early exit if derivs is large enough
            unsigned int old_df { this->order_ };
            if (this->degree_ < derivs) {
                if (complete_basis) {
                    return arma::zeros(this->x_.n_elem, old_df);
                }
                if (old_df == 1) {
                    throw std::range_error("No column left in the matrix.");
                }
                return arma::zeros(this->x_.n_elem, old_df - 1);
            }
            // back up current results if necessary
            bool backup_basis { this->is_basis_latest_ };
            rmat old_basis;
            if (backup_basis) {
                old_basis = this->poly_basis_;
            }
            // get basis matrix for (degree - derivs)
            this->set_degree(this->degree_ - derivs);
            rmat d_mat { this->basis(true) };
            // restore
            this->set_degree(this->degree_ + derivs);
            this->is_basis_latest_ = backup_basis;
            if (backup_basis) {
                this->poly_basis_ = old_basis;
            }
            // add zero columns
            d_mat = add_zero_cols(d_mat, old_df - d_mat.n_cols);
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
        inline virtual rmat integral(const bool complete_basis = true)
        {
            // back up current results if necessary
            bool backup_basis { this->is_basis_latest_ };
            rmat old_basis;
            if (backup_basis) {
                old_basis = this->poly_basis_;
            }
            // get basis matrix for (degree + 1) with intercept
            this->set_degree(this->order_);
            rmat i_mat { this->basis(false) };
            // restore
            this->set_degree(this->degree_ - 1);
            this->is_basis_latest_ = backup_basis;
            if (backup_basis) {
                this->poly_basis_ = old_basis;
            }
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
