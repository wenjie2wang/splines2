##
## R package splines2 by Wenjie Wang and Jun Yan
## Copyright (C) 2016-2021
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

##' Evaluate a Spline Basis at specified points
##'
##' This function evaluates a predefined spline basis at a (new) given \code{x}.
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


##' @rdname predict
##' @export
predict.bSpline2 <- function(object, newx, ...)
{
    if (missing(newx))
        return(object)
    ## checks if key attributes still exist
    check_attr(object, c("x", "degree", "knots", "Boundary.knots", "intercept"))
    do.call(bSpline, c(list(x = newx), pred_attr(object)))
}


##' @rdname predict
##' @export
predict.dbs <- function(object, newx, ...)
{
    if (missing(newx))
        return(object)
    ## checks if key attributes still exist
    ## checks if key attributes still exist
    check_attr(object, c("x", "degree", "derivs",
                         "knots", "Boundary.knots", "intercept"))
    do.call(dbs, c(list(x = newx), pred_attr(object)))
}


##' @rdname predict
##' @export
predict.ibs <- function(object, newx, ...)
{
    if (missing(newx))
        return(object)
    ## checks if key attributes still exist
    check_attr(object, c("x", "degree", "knots", "Boundary.knots", "intercept"))
    do.call(ibs, c(list(x = newx), pred_attr(object)))
}


##' @rdname predict
##' @export
predict.mSpline <- function(object, newx, ...)
{
    if (missing(newx))
        return(object)
    ## checks if key attributes still exist
    check_attr(object, c("x", "degree", "derivs", "integral", "periodic",
                         "knots", "Boundary.knots", "intercept"))
    do.call(mSpline, c(list(x = newx), pred_attr(object)))
}


##' @rdname predict
##' @export
predict.iSpline <- function(object, newx, ...)
{
    if (missing(newx))
        return(object)
    ## checks if key attributes still exist
    check_attr(object, c("x", "degree", "derivs",
                         "knots", "Boundary.knots", "intercept"))
    do.call(iSpline, c(list(x = newx), pred_attr(object)))
}


##' @rdname predict
##' @export
predict.cSpline <- function(object, newx, ...)
{
    if (missing(newx))
        return(object)
    ## checks if key attributes still exist
    check_attr(object, c("x", "degree", "derivs", "scale",
                         "knots", "Boundary.knots", "intercept"))
    do.call(cSpline, c(list(x = newx), pred_attr(object)))
}


##' @rdname predict
##' @export
predict.bernsteinPoly <- function(object, newx, ...)
{
    if (missing(newx))
        return(object)
    ## checks if key attributes still exist
    check_attr(object, c("x", "degree", "derivs", "integral",
                         "Boundary.knots", "intercept"))
    do.call(bernsteinPoly, c(list(x = newx), pred_attr(object)))
}


##' @rdname predict
##' @export
predict.naturalSpline <- function(object, newx, ...)
{
    if (missing(newx))
        return(object)
    ## checks if key attributes still exist
    check_attr(object, c("x", "derivs", "integral",
                         "knots", "Boundary.knots", "intercept"))
    do.call(naturalSpline, c(list(x = newx), pred_attr(object)))
}
