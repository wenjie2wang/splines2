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

#ifndef SPLINES2_ISPLINE_H
#define SPLINES2_ISPLINE_H

#include <stdexcept>

#include "common.h"
#include "utils.h"
#include "SplineBase.h"
#include "BSpline.h"
#include "MSpline.h"

namespace splines2 {

    // define a class for M-splines
    class ISpline : public SplineBase
    {
        // inherits constructors
        using SplineBase::SplineBase;

    public:
        // function members

        //! Compute I-spline basis
        //!
        //! @param complete_basis A `bool` value indicating whether to return a
        //! complete spline basis
        //!
        //! @return arma::mat
        inline virtual rmat basis(const bool complete_basis = true)
        {
            MSpline msp_obj { this };
            rmat out { msp_obj.integral(true) };
            if (complete_basis) {
                return out;
            }
            return mat_wo_col1(out);
        }

        inline virtual rmat derivative(const unsigned int derivs = 1,
                                       const bool complete_basis = true)
        {
            if (derivs == 0) {
                throw std::range_error(
                    "'derivs' has to be a positive integer.");
            }
            MSpline msp_obj { this };
            if (derivs == 1) {
                return msp_obj.basis(complete_basis);
            }
            return msp_obj.derivative(derivs - 1, complete_basis);
        }

        inline virtual rmat integral(const bool complete_basis = true)
        {
            BSpline bsp_obj { this };
            bsp_obj.set_degree(degree_ + 1);
            rmat i_mat { bsp_obj.integral(false) };
            update_x_index();
            // for each row of i_mat
            for (size_t i {0}; i < x_.n_elem; ++i) {
                size_t k2 { x_index_(i) + degree_ };
                arma::rowvec numer { i_mat(i, arma::span(0, k2)) };
                numer = rev_cum_sum(numer);
                for (size_t j {0}; j < i_mat.n_cols; ++j) {
                    if (j > k2) {
                        i_mat(i, j) = 0.0;
                    } else {
                        i_mat(i, j) = numer(j);
                    }
                }
            }
            // remove the first column if needed
            if (complete_basis) {
                return i_mat;
            }
            return mat_wo_col1(i_mat);
        }

    };
}  // splines2


#endif /* SPLINES2_ISPLINE_H */
