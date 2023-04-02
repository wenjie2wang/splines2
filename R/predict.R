##
## R package splines2 by Wenjie Wang and Jun Yan
## Copyright (C) 2016-2023
##
## This file is part of the R package splines2.
##
## The R package splines2 is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package splines2 is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##

##' Evaluate Spline Basis Functions at Specified Points
##'
##' This function evaluates the given spline basis functions at the specified
##' \code{x}.
##'
##' These are methods for the generic function \code{predict} for objects
##' inheriting from class \code{bSpline2}, \code{ibs}, \code{mSpline},
##' \code{iSpline}, \code{cSpline}, \code{naturalSpline}, or
##' \code{bernsteinPoly}.  If \code{newx} is not given, the function returns the
##' input object.
##'
##' @name predict
##' @param object Objects of class \code{bSpline2}, \code{ibs}, \code{mSpline},
##'     \code{iSpline}, \code{cSpline}, \code{bernsteinPoly} or
##'     \code{naturalSpline} with attributes describing \code{knots},
##'     \code{degree}, etc.
##' @param newx The \code{x} values at which evaluations are required.
##' @param ... Optional arguments that are not used.
##'
##' @return
##' An object just like the \code{object} input, except evaluated at
##' the new values of \code{x}.
##'
##' @example inst/examples/ex-predict.R
##'
##' @importFrom stats predict
NULL


## the helper function for predict
helper_predict <- function(object, newx, fun, key_attr)
{
    if (missing(newx))
        return(object)
    ## checks if key attributes still exist
    check_attr(object, key_attr)
    do.call(fun, c(list(x = newx), pred_attr(object)))
}


##' @rdname predict
##' @export
predict.BSpline <- function(object, newx, ...)
{
    helper_predict(
        object,
        newx,
        bSpline,
        key_attr = c("x", "degree", "knots", "Boundary.knots",
                     "intercept", "periodic", "derivs", "integral")
    )
}

##' @rdname predict
##' @export
predict.MSpline <- function(object, newx, ...)
{
    helper_predict(
        object,
        newx,
        mSpline,
        key_attr = c("x", "degree", "knots", "Boundary.knots",
                     "intercept", "periodic", "derivs", "integral")
    )
}


##' @rdname predict
##' @export
predict.ISpline <- function(object, newx, ...)
{
    helper_predict(
        object,
        newx,
        iSpline,
        key_attr = c("x", "degree", "knots", "Boundary.knots",
                     "intercept", "derivs")
    )
}


##' @rdname predict
##' @export
predict.CSpline <- function(object, newx, ...)
{
    helper_predict(
        object,
        newx,
        cSpline,
        key_attr = c("x", "degree", "knots", "Boundary.knots",
                     "intercept", "derivs", "scale")
    )
}


##' @rdname predict
##' @export
predict.BernsteinPoly <- function(object, newx, ...)
{
    helper_predict(
        object,
        newx,
        bernsteinPoly,
        key_attr = c("x", "degree", "Boundary.knots",
                     "intercept", "derivs", "integral")
    )
}


##' @rdname predict
##' @export
predict.NaturalSpline <- function(object, newx, ...)
{
    helper_predict(
        object,
        newx,
        naturalSpline,
        key_attr = c("x", "knots", "Boundary.knots",
                     "intercept", "derivs", "integral")
    )
}
