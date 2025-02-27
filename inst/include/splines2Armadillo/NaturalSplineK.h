//
// R package splines2 by Wenjie Wang and Jun Yan
// Copyright (C) 2016-2025
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

#ifndef SPLINES2_NATURALSPLINEK_H
#define SPLINES2_NATURALSPLINEK_H

#include <stdexcept>

#include "common.h"
#include "utils.h"
#include "NaturalSpline.h"

namespace splines2 {

    // define a class for Natural cubic splines
    class NaturalSplineK : public NaturalSpline
    {
    protected:

        // compute the transformation matrix so that only one of the basis
        // functions attains one at the internal and boundary knots
        inline void set_null_colvecs(const bool standardize = true) override
        {
            // get natural spline basis functions at knots
            NaturalSpline nsp_obj { this };
            rvec flat_knots {
                get_simple_knot_sequence(nsp_obj.get_internal_knots(),
                                         nsp_obj.get_boundary_knots(), 1)
            };
            nsp_obj.set_x(flat_knots);
            rmat kmat { nsp_obj.basis(true) };
            null_colvecs_ = nsp_obj.get_transform_matrix() * arma::inv(kmat);
        }

    public:
        // the default constructor
        NaturalSplineK() {}

        using NaturalSpline::NaturalSpline;

    };


}  // splines2


#endif /* SPLINES2_NATURALSPLINEK_H */
