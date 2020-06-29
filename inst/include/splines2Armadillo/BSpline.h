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

#ifndef SPLINES2_BSPLINE_H
#define SPLINES2_BSPLINE_H

#include <stdexcept>
#include <memory>

#include "aliases.h"
#include "utils.h"
#include "SplineBase.h"

namespace splines2arma {

    // define a class for B-splines
    class BSpline : public SplineBase
    {
        // inherits constructors
        using SplineBase::SplineBase;

    public:
        // function members

        //! Compute B-spline basis
        //!
        //! @param complete_basis A `bool` value indicating whether to return a
        //! complete spline basis
        //!
        //! @return arma::mat
        inline rmat basis(const bool complete_basis = true)
        {
            // early exit if latest
            if (this->is_basis_latest_) {
                if (complete_basis) {
                    return this->spline_basis_;
                }
                // else
                return mat_wo_col1(this->spline_basis_);
            }
            // else do the generation
            this->update_x_index();
            // define output matrix
            rmat b_mat {
                arma::zeros(this->x_.n_elem, this->spline_df_)
            };
            // generate bases of degree 0
            for (size_t i {0}; i < x_.n_elem; ++i) {
                b_mat(i, this->x_index_(i)) = 1;
            }
            // only need knot sequence for degree > 0
            if (degree_ > 0) {
                this->update_knot_sequence();
            }
            // main loop
            for (unsigned int k {1}; k <= degree_; ++k) {
                const unsigned int k_offset { degree_ - k };
                // use the Cox-de Boor recursive formula
                for (size_t i {0}; i < x_.n_elem; ++ i) {
                    double saved { 0 };
                    // for each x, at most "order" columns are not zero
                    // basis(x) is not zero from t_ii to t_{ii+k+1}
                    // where ii is index of x in terms of bases
                    for (size_t j {0}; j < k; ++j) {
                        size_t j_index { x_index_(i) + j };
                        size_t i1 { j_index + k_offset + 1 };
                        size_t i2 { j_index + order_ };
                        double den { knot_sequence_(i2) - knot_sequence_(i1) };
                        double term { b_mat(i, j_index) / den };
                        b_mat(i, j_index) = saved +
                            (knot_sequence_(i2) - x_(i)) * term;
                        saved = (x_(i) - knot_sequence_(i1)) * term;
                    }
                    b_mat(i, x_index_(i) + k) = saved;
                }
            }
            // about to return
            this->spline_basis_ = b_mat;
            this->is_basis_latest_ = true;
            if (complete_basis) {
                return b_mat;
            }
            // else
            return mat_wo_col1(b_mat);
        }

        // derivatives of B-splines
        inline rmat deriv(
            const unsigned int derivs = 1,
            const bool complete_basis = true
            )
        {
            if (derivs == 0) {
                throw std::range_error(
                    "'derivs' has to be a positive integer."
                    );
            }
            // early exit if derivs is large enough
            unsigned int old_df { this->spline_df_ };
            if (this->degree_ < derivs) {
                if (complete_basis) {
                    return arma::zeros(this->x_.n_elem, old_df);
                }
                return arma::zeros(this->x_.n_elem, old_df - 1);
            }
            // back up current results if necessary
            bool backup_basis { this->is_basis_latest_ };
            bool backup_knot_sequence { this->is_knot_sequence_latest_ };
            rmat old_basis;
            rvec old_knot_sequence;
            if (backup_basis) {
                old_basis = this->spline_basis_;
            }
            if (backup_knot_sequence) {
                old_knot_sequence = this->knot_sequence_;
            }
            // get basis matrix for (degree - derivs)
            this->set_degree(this->degree_ - derivs);
            rmat d_mat { this->basis(true) };
            // restore
            this->set_degree(this->degree_ + derivs);
            if (backup_basis) {
                this->is_basis_latest_ = true;
                this->spline_basis_ = old_basis;
            }
            if (backup_knot_sequence) {
                this->is_knot_sequence_latest_ = true;
                this->knot_sequence_ = old_knot_sequence;
            } else {
                this->update_knot_sequence();
            }
            // add zero columns
            d_mat = add_zero_cols(d_mat, old_df - d_mat.n_cols);
            // main loop
            for (unsigned int k {1}; k <= derivs; ++k) {
                const unsigned int k_offset { derivs - k };
                const size_t numer { degree_ - k_offset };
                for (size_t i {0}; i < x_.n_elem; ++ i) {
                    double saved { 0 };
                    for (size_t j {0}; j < numer; ++j) {
                        size_t j_index { x_index_(i) + j };
                        size_t i1 { j_index + k_offset + 1 };
                        size_t i2 { j_index + order_ };
                        double den { knot_sequence_(i2) - knot_sequence_(i1) };
                        double term { numer * d_mat(i, j_index) / den };
                        d_mat(i, j_index) = saved - term;
                        saved = term;
                    }
                    d_mat(i, x_index_(i) + numer) = saved;
                }
            }
            return d_mat;
        }

    };

}

#endif
