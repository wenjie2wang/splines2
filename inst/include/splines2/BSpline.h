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

#ifndef SPLINE2_B_SPLINE_H
#define SPLINE2_B_SPLINE_H

#include <stdexcept>

#include "common.h"
#include "utils.h"
#include "SplineBase.h"

namespace splines2 {

    // define a class for B-splines
    class BSpline : public SplineBase
    {
        using SplineBase::SplineBase;
    public:
        // function members
        rmat basis(const rvec& new_x = rvec(),
                   const bool update_x = true,
                   const bool complete_basis = true)
        {
            // define some variables
            bool is_x_new { new_x.size() > 0 };
            // early exit if no new x and latest
            if (! is_x_new && this->is_basis_latest) {
                if (complete_basis) {
                    return this->full_basis;
                }
                // else
                return mat_wo_col1(this->full_basis);
            }
            // else do the generation
            rvec old_x { rvec() };
            ivec old_x_index { ivec() };
            // if a meaningful nex_x is specified
            if (is_x_new) {
                if (! update_x) {
                    // in case the x_index is not latest
                    this->update_x_index();
                    // save the previous x and x_index
                    old_x = this->x;
                    old_x_index = this->x_index;
                }
                // else, update x
                this->x = new_x;
            }
            this->update_x_index();
            // define output matrix
            rmat b_mat {
                rmat(this->x.size(), this->internal_knots.size() + this->order)
            };
            // generate bases of degree 0
            for (long i {0}; i < x.size(); ++i) {
                b_mat(i, x_index(i)) = 1;
            }
            // only need knot sequence for degree > 0
            if (degree > 0) {
                this->update_knot_sequence();
            }
            // main loop
            for (unsigned int k {1}; k <= degree; ++k) {
                const unsigned int k_offset { degree - k };
                // use the Cox-de Boor recursive formula
                for (long i {0}; i < x.size(); ++ i) {
                    double saved { 0 };
                    // for each x, at most "order" columns are not zero
                    // x is not zero from t_ii to t_{ii+k+1}
                    // where ii is index of basis
                    for (size_t j {0}; j < k; ++j) {
                        unsigned long j_index { x_index(i) + j };
                        unsigned long i1 { j_index + k_offset + 1 };
                        unsigned long i2 { j_index + order };
                        double den { knot_sequence(i2) - knot_sequence(i1) };
                        double term { b_mat(i, j_index) / den };
                        b_mat(i, j_index) = saved +
                            (knot_sequence(i2) - x(i)) * term;
                        saved = (x(i) - knot_sequence(i1)) * term;
                    }
                    b_mat(i, x_index(i) + k) = saved;
                }
            }
            // about to return
            if (! is_x_new) {
                this->is_basis_latest = true;
            } else if (! update_x) {
                // set old x and index back for new_x
                this->x = old_x;
                this->x_index = old_x_index;
            }
            if (complete_basis) {
                return b_mat;
            }
            // else
            return mat_wo_col1(b_mat);
        }

        // return bases for pre-set x
        rmat basis(const bool complete_basis = true)
        {
            return basis(rvec(), false, complete_basis);
        }



    };

}




#endif
