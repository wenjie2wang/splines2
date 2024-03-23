##
## R package splines2 by Wenjie Wang and Jun Yan
## Copyright (C) 2016-2024
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

##' Derivatives of Spline Basis Functions
##'
##' Returns derivatives of given order for the given spline basis functions.
##'
##' At knots, the derivative is defined to be the right derivative except at the
##' right boundary knot. By default, the function returns the first derivatives.
##' For derivatives of order greater than one, nested function calls such as
##' \code{deriv(deriv(expr))} are supported but not recommended.  For a better
##' performance, argument \code{derivs} should be specified instead.
##'
##' This function is designed for objects produced by this package.  It
##' internally extracts necessary specification about the spline/polynomial
##' basis matrix from its attributes. Therefore, the function will not work if
##' the key attributes are not available after some operations.
##'
##' @name deriv
##'
##' @param expr Objects of class \code{bSpline2}, \code{ibs}, \code{mSpline},
##'     \code{iSpline}, \code{cSpline}, \code{bernsteinPoly} or
##'     \code{naturalSpline} with attributes describing \code{knots},
##'     \code{degree}, etc.
##' @param derivs A positive integer specifying the order of derivatives. By
##'     default, it is \code{1L} for the first derivatives.
##' @param ... Optional arguments that are not used.
##'
##' @return A numeric matrix of the same dimension with the input \code{expr}.
##'
##' @example inst/examples/ex-deriv.R
##'
##' @importFrom stats deriv
NULL


##' @rdname deriv
##' @export
deriv.BSpline <- function(expr, derivs = 1L, ...)
{
    ## quick check on derivs
    derivs <- as.integer(derivs)
    if (derivs < 1L) {
        stop("The 'derivs' has to be a positive integer.")
    }
    ## checks if key attributes still exist
    check_attr(expr, c("x", "degree", "knots", "Boundary.knots", "intercept",
                       "periodic", "derivs", "integral"))
    attr(expr, "derivs") <- attr(expr, "derivs") + derivs
    do.call(bSpline, attributes(expr))
}


##' @rdname deriv
##' @export
deriv.MSpline <- function(expr, derivs = 1L, ...)
{
    ## quick check on derivs
    derivs <- as.integer(derivs)
    if (derivs < 1L) {
        stop("The 'derivs' has to be a positive integer.")
    }
    ## checks if key attributes still exist
    check_attr(expr, c("x", "degree", "knots", "Boundary.knots", "intercept",
                       "periodic", "derivs", "integral"))
    attr(expr, "derivs") <- attr(expr, "derivs") + derivs
    do.call(mSpline, attributes(expr))
}


##' @rdname deriv
##' @export
deriv.ISpline <- function(expr, derivs = 1L, ...)
{
    ## quick check on derivs
    derivs <- as.integer(derivs)
    if (derivs < 1L) {
        stop("The 'derivs' has to be a positive integer.")
    }
    ## checks if key attributes still exist
    check_attr(expr, c("x", "degree", "derivs",
                       "knots", "Boundary.knots", "intercept"))
    attr(expr, "derivs") <- derivs - 1L
    do.call(mSpline, attributes(expr))
}


##' @rdname deriv
##' @export
deriv.CSpline <- function(expr, derivs = 1L, ...)
{
    ## quick check on derivs
    derivs <- as.integer(derivs)
    if (derivs < 1L) {
        stop("The 'derivs' has to be a positive integer.")
    }
    ## checks if key attributes still exist
    check_attr(expr, c("x", "degree", "derivs",
                       "knots", "Boundary.knots", "intercept"))
    scl <- attr(expr, "scale")
    if (scl) {
        ## if scaled
        attr(expr, "derivs") <- attr(expr, "derivs") + derivs
        do.call(cSpline, attributes(expr))
    } else {
        ## if not scaled (then "derivs" must be 0 in cSpline call)
        derivs <- as.integer(derivs)
        if (derivs == 1L) {
            return(do.call(iSpline, attributes(expr)))
        }
        if (derivs == 2L) {
            return(do.call(mSpline, attributes(expr)))
        }
        attr(expr, "derivs") <- derivs - 2L
        do.call(mSpline, attributes(expr))
    }
}


##' @rdname deriv
##' @export
deriv.BernsteinPoly <- function(expr, derivs = 1L, ...)
{
    ## quick check on derivs
    derivs <- as.integer(derivs)
    if (derivs < 1L) {
        stop("The 'derivs' has to be a positive integer.")
    }
    ## checks if key attributes still exist
    check_attr(expr, c("x", "degree", "derivs", "integral",
                       "Boundary.knots", "intercept"))
    attr(expr, "derivs") <- attr(expr, "derivs") + derivs
    do.call(bernsteinPoly, attributes(expr))
}


##' @rdname deriv
##' @export
deriv.NaturalSpline <- function(expr, derivs = 1L, ...)
{
    ## quick check on derivs
    derivs <- as.integer(derivs)
    if (derivs < 1L) {
        stop("The 'derivs' has to be a positive integer.")
    }
    ## checks if key attributes still exist
    check_attr(expr, c("x", "derivs", "integral", "knots",
                       "Boundary.knots", "trim",  "intercept"))
    attr(expr, "derivs") <- attr(expr, "derivs") + derivs
    do.call(naturalSpline, attributes(expr))
}

##' @rdname deriv
##' @export
deriv.NaturalSplineK <- function(expr, derivs = 1L, ...)
{
    ## quick check on derivs
    derivs <- as.integer(derivs)
    if (derivs < 1L) {
        stop("The 'derivs' has to be a positive integer.")
    }
    ## checks if key attributes still exist
    check_attr(expr, c("x", "derivs", "integral", "knots",
                       "Boundary.knots", "trim",  "intercept"))
    attr(expr, "derivs") <- attr(expr, "derivs") + derivs
    do.call(nsk, attributes(expr))
}
