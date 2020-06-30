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

#ifndef SPLINE2_ISPLINE_H
#define SPLINE2_ISPLINE_H

#include <memory>
#include <stdexcept>

#include "aliases.h"
#include "utils.h"
#include "SplineBase.h"
#include "MSpline.h"

namespace splines2arma {

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
            // else do generation
            MSpline msp_obj { MSpline(this) };
            this->spline_basis_ = msp_obj.integral(complete_basis);
            this->is_basis_latest_ = true;
            return this->spline_basis_;
        }

        inline rmat derivative(const unsigned int derivs = 1,
                               const bool complete_basis = true)
        {
            MSpline msp_obj { MSpline(this) };
            if (derivs == 1) {
                return msp_obj.basis(complete_basis);
            }
            return msp_obj.derivative(derivs - 1, complete_basis);
        }

        inline rmat integral(const bool complete_basis = true)
        {
            return rmat();
        }

    };
}  // splines2arma


#endif /* SPLINE2_ISPLINE_H */
