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

        // indices of x placed outside of boundary (left/right)
        bool is_x_outside_latest_ = false;
        uvec x_outside_left_;
        uvec x_outside_right_;

        using SplineBase::set_degree;
        using SplineBase::set_order;
        using SplineBase::SplineBase;

        // get null space vector for the second derivatives
        // of B-spline basis on boundary knotsn
        // the results depend on knot sequence only
        inline void set_null_colvecs(bool standardize = true)
        {
            null_colvecs_ = arma::zeros<arma::mat>(spline_df_ + 2, spline_df_);
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
                    null_colvecs_(spline_df_ - i + 1, 1) = 1.0;
                }
                null_colvecs_(1, 2) = 1.0;
                null_colvecs_(2, 2) = 1.0 +
                    (internal_knots_(1) - boundary_knots_(0)) /
                    (internal_knots_(0) - boundary_knots_(0));
                null_colvecs_(spline_df_ - 1, 3) = 1.0 +
                    (boundary_knots_(1) - internal_knots_(n_knots - 2)) /
                    (boundary_knots_(1) - internal_knots_(n_knots - 1));
                null_colvecs_(spline_df_, 3) = 1.0;
                if (spline_df_ > 4) {
                    for (size_t j {3}; j < spline_df_ - 1; ++j) {
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

        // update x index for outside
        inline void update_x_outside()
        {
            if (! is_x_outside_latest_) {
                x_outside_left_ = arma::find(x_ < boundary_knots_(0));
                x_outside_right_ = arma::find(x_ > boundary_knots_(1));
                is_x_outside_latest_ = true;
            }
        }

        // compute spline df
        inline void update_spline_df()
        {
            spline_df_ = internal_knots_.n_elem + 2;
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
            update_spline_df();
            this->set_null_colvecs();
            update_x_outside();
        }

        // constructor with specificied internal_knots
        NaturalSpline(const rvec& x,
                      const rvec& internal_knots,
                      const rvec& boundary_knots = rvec())
        {
            x_ = x;
            degree_ = 3;
            order_ = 4;
            simplify_knots(internal_knots, boundary_knots);
            update_spline_df();
            update_x_outside();
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
            unsigned int n_internal_knots { spline_df_ - 2 };
            if (n_internal_knots == 0) {
                simplify_knots(rvec(), boundary_knots);
            } else {
                rvec prob_vec { arma::linspace(0, 1, n_internal_knots + 2) };
                prob_vec = prob_vec.subvec(1, n_internal_knots);
                simplify_knots(rvec(), boundary_knots);
                // get quantiles of x within boundary only
                rvec x_inside { get_inside_x(x, boundary_knots_) };
                rvec internal_knots { arma_quantile(x_inside, prob_vec) };
                simplify_knots(internal_knots);
            }
            update_x_outside();
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
            this->set_null_colvecs();
            BSpline bs_obj { this };
            rmat bsMat { bs_obj.basis(true) };
            // precess x outside of boundary
            update_x_outside();
            if (x_outside_left_.n_elem > 0 || x_outside_right_.n_elem > 0) {
                BSpline bs_tmp;
                bs_tmp.set_degree(3);
                bs_tmp.set_internal_knots(internal_knots_);
                bs_tmp.set_boundary_knots(boundary_knots_);
                if (x_outside_left_.n_elem > 0) {
                    bs_tmp.set_x(boundary_knots_(0));
                    rmat tt1 { bs_tmp.basis(true) };
                    rmat tt2 { bs_tmp.derivative(true) };
                    for (size_t k {0}; k < x_outside_left_.n_elem; ++k) {
                        size_t idx { x_outside_left_(k) };
                        bsMat.row(idx) = tt1 +
                            tt2 * (x_(idx) - boundary_knots_(0));
                    }
                }
                if (x_outside_right_.n_elem > 0) {
                    bs_tmp.set_x(boundary_knots_(1));
                    rmat tt1 { bs_tmp.basis(true) };
                    rmat tt2 { bs_tmp.derivative(true) };
                    for (size_t k {0}; k < x_outside_right_.n_elem; ++k) {
                        size_t idx { x_outside_right_(k) };
                        bsMat.row(idx) = tt1 +
                            tt2 * (x_(idx) - boundary_knots_(1));
                    }
                }
            }
            // apply null space
            bsMat *= null_colvecs_;
            if (complete_basis) {
                return bsMat;
            }
            // else
            return mat_wo_col1(bsMat);
        }

        inline rmat derivative(const unsigned int derivs = 1,
                               const bool complete_basis = true)
        {
            this->set_null_colvecs();
            BSpline bs_obj { this };
            rmat bsMat { bs_obj.derivative(derivs, true) };
            // precess x outside of boundary
            update_x_outside();
            if (x_outside_left_.n_elem > 0 || x_outside_right_.n_elem > 0) {
                if (derivs > 1) {
                    arma::rowvec zero_row {
                        arma::zeros<arma::rowvec>(bsMat.n_cols)
                    };
                    if (x_outside_left_.n_elem > 0) {
                        for (size_t k {0}; k < x_outside_left_.n_elem; ++k) {
                            size_t idx { x_outside_left_(k) };
                            bsMat.row(idx) = zero_row;
                        }
                    }
                    if (x_outside_right_.n_elem > 0) {
                        for (size_t k {0}; k < x_outside_right_.n_elem; ++k) {
                            size_t idx { x_outside_right_(k) };
                            bsMat.row(idx) = zero_row;
                        }
                    }
                } else {
                    // derivs = 1
                    BSpline bs_tmp;
                    bs_tmp.set_degree(3);
                    bs_tmp.set_internal_knots(internal_knots_);
                    bs_tmp.set_boundary_knots(boundary_knots_);
                    if (x_outside_left_.n_elem > 0) {
                        bs_tmp.set_x(boundary_knots_(0));
                        rmat tt2 { bs_tmp.derivative(true) };
                        for (size_t k {0}; k < x_outside_left_.n_elem; ++k) {
                            size_t idx { x_outside_left_(k) };
                            bsMat.row(idx) = tt2;
                        }
                    }
                    if (x_outside_right_.n_elem > 0) {
                        bs_tmp.set_x(boundary_knots_(1));
                        rmat tt2 { bs_tmp.derivative(true) };
                        for (size_t k {0}; k < x_outside_right_.n_elem; ++k) {
                            size_t idx { x_outside_right_(k) };
                            bsMat.row(idx) = tt2;
                        }
                    }
                }
            }
            // apply null space
            rmat out { bsMat * null_colvecs_ };
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
            rmat bsMat { bs_obj.integral(true) };
            // precess x outside of boundary
            update_x_outside();
            if (x_outside_left_.n_elem > 0 || x_outside_right_.n_elem > 0) {
                // integrate from left boundary
                if (x_outside_left_.n_elem > 0) {
                    arma::rowvec zero_row {
                        arma::zeros<arma::rowvec>(bsMat.n_cols)
                    };
                    for (size_t k {0}; k < x_outside_left_.n_elem; ++k) {
                        size_t idx { x_outside_left_(k) };
                        bsMat.row(idx) = zero_row;
                    }
                }
                if (x_outside_right_.n_elem > 0) {
                    BSpline bs_tmp;
                    bs_tmp.set_degree(3);
                    bs_tmp.set_internal_knots(internal_knots_);
                    bs_tmp.set_boundary_knots(boundary_knots_);
                    bs_tmp.set_x(boundary_knots_(1));
                    arma::rowvec right_row {
                        bs_obj.set_x(boundary_knots_(1))->integral(true).row(0)
                    };
                    rmat tt1 { bs_tmp.basis(true) };
                    rmat tt2 { bs_tmp.derivative(true) };
                    for (size_t k {0}; k < x_outside_right_.n_elem; ++k) {
                        size_t idx { x_outside_right_(k) };
                        double tmp { x_(idx) - boundary_knots_(1) };
                        bsMat.row(idx) = right_row + tt1 * tmp +
                            0.5 * tt2 * tmp * tmp;
                    }
                }
            }
            // apply null space
            rmat out { bsMat * null_colvecs_ };
            if (complete_basis) {
                return out;
            }
            // else
            return mat_wo_col1(out);
        }

        // re-define some "setter" functions
        inline virtual NaturalSpline* set_x(const rvec& x)
        {
            x_ = x;
            is_x_index_latest_ = false;
            is_x_outside_latest_ = false;
            return this;
        }
        inline virtual NaturalSpline* set_x(const double x)
        {
            x_ = num2vec(x);
            is_x_index_latest_ = false;
            is_x_outside_latest_ = false;
            return this;
        }
        inline virtual NaturalSpline* set_boundary_knots(
            const rvec& boundary_knots
            )
        {
            simplify_knots(internal_knots_, boundary_knots);
            is_knot_sequence_latest_ = false;
            is_x_index_latest_ = false;
            is_x_outside_latest_ = false;
            return this;
        }


    };


}  // splines2


#endif /* SPLINES2_NATURALSPLINE_H */
