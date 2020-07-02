##
## R package splines2 by Wenjie Wang and Jun Yan
## Copyright (C) 2016-2020
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

##' Evaluate a Spline Basis
##'
##' This function evaluates a predefined spline basis at (new) given values.
##'
##' These are methods for the generic function \code{predict} for objects
##' inheriting from class \code{bSpline2}, \code{ibs}, \code{mSpline},
##' \code{iSpline}, or \code{cSpline}.  If \code{newx} is not given, the
##' function returns the input object.  For object returned by function
##' \code{\link{cSpline}}, the \code{mSpline} and \code{iSpline} objects shipped
##' in attributes should not be evaluated by this function if \code{rescale} is
##' \code{TRUE}.  See \code{\link{cSpline}} for details.
##'
##' @name predict
##' @param object Objects of class \code{bSpline2}, \code{ibs}, \code{mSpline},
##'     \code{iSpline}, or \code{cSpline} having attributes describing
##'     \code{knots}, \code{degree}, etc.
##' @param newx The \code{x} values at which evaluations are required.
##' @param ... Optional argument for future usage.
##'
##' @return
##' An object just like the \code{object} input, except evaluated at
##' the new values of \code{x}.
##'
##' @example inst/examples/ex-predict.R
##'
##' @seealso
##' \code{\link{bSpline}} for B-splines;
##' \code{\link{ibs}} for integral of B-splines;
##' \code{\link{dbs}} for derivative of B-splines;
##' \code{\link{mSpline}} for M-splines;
##' \code{\link{iSpline}} for I-splines;
##' \code{\link{cSpline}} for C-splines.
##'
##' @importFrom stats predict
NULL


##' @rdname predict
##' @export
predict.bSpline2 <- function(object, newx, ...)
{
    if (missing(newx))
        return(object)
    a <- c(list(x = newx),
           attributes(object)[c("degree", "knots", "Boundary.knots",
                                "intercept")])
    do.call("bSpline", a)
}


##' @rdname predict
##' @export
predict.ibs <- function(object, newx, ...)
{
    if (missing(newx))
        return(object)
    a <- c(list(x = newx),
           attributes(object)[c("degree", "knots", "Boundary.knots",
                                "intercept")])
    do.call("ibs", a)
}


##' @rdname predict
##' @export
predict.dbs <- function(object, newx, ...)
{
    if (missing(newx))
        return(object)
    a <- c(list(x = newx),
           attributes(object)[c("degree", "knots", "Boundary.knots",
                                "intercept", "derivs")])
    do.call("dbs", a)
}


##' @rdname predict
##' @export
predict.mSpline <- function(object, newx, ...)
{
    if (missing(newx))
        return(object)
    a <- c(list(x = newx),
           attributes(object)[c("degree", "knots", "Boundary.knots",
                                "intercept", "derivs")])
    do.call("mSpline", a)
}


##' @rdname predict
##' @export
predict.iSpline <- function(object, newx, ...)
{
    if (missing(newx))
        return(object)
    a <- c(list(x = newx),
           attributes(object)[c("degree", "knots", "Boundary.knots",
                                "intercept", "derivs")])
    do.call("iSpline", a)
}


##' @rdname predict
##' @export
predict.cSpline <- function(object, newx, ...)
{
    if (missing(newx))
        return(object)
    a <- c(list(x = newx),
           attributes(object)[c("degree", "knots", "Boundary.knots",
                                "intercept", "scale")])
    do.call("cSpline", a)
}
