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

#ifndef SPLINES2_NATURALSPLINE_H
#define SPLINES2_NATURALSPLINE_H

#include <stdexcept>

#include "common.h"
#include "utils.h"
#include "SplineBase.h"
#include "BSpline.h"

namespace splines2 {

    // define a class for Natural splines
    class NaturalSpline : public SplineBase
    {
    protected:
        rmat null_colvecs_;

        using SplineBase::set_degree;
        using SplineBase::set_order;

        // get null space vector for the second derivatives
        // of B-spline basis on boundary knotsn
        // the results depend on knot sequence only
        inline void set_null_colvecs(bool standardize = true)
        {
            null_colvecs_ = arma::zeros<arma::mat>(spline_df_, spline_df_ - 2);
            size_t n_knots { internal_knots_.n_elem };
            if (n_knots == 0) {
                // without any internal knot
                null_colvecs_(0, 0) = 3.0;
                null_colvecs_(1, 0) = 2.0;
                null_colvecs_(2, 0) = 1.0;
                null_colvecs_(1, 1) = 1.0;
                null_colvecs_(2, 1) = 2.0;
                null_colvecs_(3, 1) = 3.0;
            } else if (n_knots == 1) {
                // with only one internal knot
                null_colvecs_(0, 0) = 1.0 +
                    (internal_knots_(0) - boundary_knots_(0)) /
                    (boundary_knots_(1) - boundary_knots_(0));
                null_colvecs_(1, 0) = 1.0;
                null_colvecs_(1, 1) = 1.0 /
                    (1.0 / (internal_knots_(0) - boundary_knots_(0)) + 1.0);
                null_colvecs_(2, 1) = 1.0;
                null_colvecs_(3, 1) = 1.0 /
                    (1.0 / (boundary_knots_(1) - internal_knots_(0)) + 1.0);
                null_colvecs_(3, 2) = 1.0;
                null_colvecs_(4, 2) = 1.0 +
                    (boundary_knots_(1) - internal_knots_(0)) /
                    (boundary_knots_(1) - boundary_knots_(0));
            } else {
                // with at least two internal knots
                for (size_t i {0}; i < 3; ++i) {
                    null_colvecs_(i, 0) = 1.0;
                    null_colvecs_(spline_df_ - i - 1, 1) = 1.0;
                }
                null_colvecs_(1, 2) = 1.0;
                null_colvecs_(2, 2) = 1.0 +
                    (internal_knots_(1) - boundary_knots_(0)) /
                    (internal_knots_(0) - boundary_knots_(0));
                null_colvecs_(spline_df_ - 3, 3) = 1.0 +
                    (boundary_knots_(1) - internal_knots_(n_knots - 2)) /
                    (boundary_knots_(1) - internal_knots_(n_knots - 1));
                null_colvecs_(spline_df_ - 2, 3) = 1.0;
                if (spline_df_ > 6) {
                    for (size_t j {3}; j < spline_df_ - 3; ++j) {
                        null_colvecs_(j, j + 1) = 1.0;
                    }
                }
            }
            // standardize coefficient for each column
            if (standardize) {
                size_t ncol_out { null_colvecs_.n_cols };
                for (size_t j {0}; j < ncol_out; ++j) {
                    null_colvecs_.col(j) /= arma::sum(null_colvecs_.col(j));
                }
            }
        }

    public:
        // the default constructor
        NaturalSpline() {}

        // explicit constructor
        explicit NaturalSpline(const SplineBase* pSplineBase) :
            SplineBase(pSplineBase)
        {
            degree_ = 3;
            order_ = 4;
            this->set_null_colvecs();
        }

        // constructor with specificied internal_knots
        NaturalSpline(const rvec& x,
                      const rvec& internal_knots,
                      const rvec& boundary_knots = rvec())
        {
            x_ = x;
            degree_ = 3;
            order_ = 4;
            clean_knots(internal_knots, boundary_knots);
            update_spline_df();
        }

        // constructor with specified df
        NaturalSpline(const rvec& x,
                      const unsigned int spline_df,
                      const rvec& boundary_knots = rvec())
        {
            x_ = x;
            degree_ = 3;
            order_ = 4;
            if (spline_df < 2) {
                // df has to be at least 2 (order - 2)
                throw std::range_error("The specified df was too small.");
            }
            spline_df_ = spline_df;
            // determine internal knots by spline_df and x
            unsigned int n_internal_knots { spline_df_ - order_ };
            if (n_internal_knots == 0) {
                clean_knots(rvec(), boundary_knots);
            } else {
                rvec prob_vec { arma::linspace(0, 1, n_internal_knots + 2) };
                prob_vec = prob_vec.subvec(1, n_internal_knots);
                clean_knots(rvec(), boundary_knots);
                // get quantiles of x within boundary only
                rvec x_inside { get_inside_x(x, boundary_knots_) };
                rvec internal_knots { arma_quantile(x_inside, prob_vec) };
                clean_knots(internal_knots);
            }
        }

        // function members

        //! Compute Natural-spline basis based on B-spline
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
            this->set_null_colvecs();
            BSpline bs_obj { this };
            spline_basis_ = bs_obj.basis(true) * null_colvecs_;
            is_basis_latest_ = true;
            if (complete_basis) {
                return spline_basis_;
            }
            // else
            return mat_wo_col1(spline_basis_);
        }

        inline rmat derivative(const unsigned int derivs = 1,
                               const bool complete_basis = true)
        {
            this->set_null_colvecs();
            BSpline bs_obj { this };
            rmat out { bs_obj.derivative(derivs, true) * null_colvecs_ };
            if (complete_basis) {
                return out;
            }
            // else
            return mat_wo_col1(out);
        }

        inline rmat integral(const bool complete_basis = true)
        {
            this->set_null_colvecs();
            BSpline bs_obj { this };
            rmat out { bs_obj.integral(true) * null_colvecs_ };
            if (complete_basis) {
                return out;
            }
            // else
            return mat_wo_col1(out);
        }

    };


}  // splines2


#endif /* SPLINES2_NATURALSPLINE_H */
