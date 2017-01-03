################################################################################
##
##   R package splines2 by Wenjie Wang and Jun Yan
##   Copyright (C) 2016-2017
##
##   This file is part of the R package splines2.
##
##   The R package splines2 is free software: You can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   any later version (at your option). See the GNU General Public License
##   at <http://www.gnu.org/licenses/> for details.
##
##   The R package splines2 is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##
################################################################################


##' Derivative of Splines
##'
##' Obtain derivative of spline bases.
##'
##' @name deriv
##'
##' @param expr \code{spline2} objects generated from this package.
##' @param derivs A positive integer indicating the number of
##' derivatives.
##' @param ... Other arguments for further usage.
##'
##' @importFrom stats deriv
NULL


##' @rdname deriv
##' @export
deriv.ibs <- function(expr, derivs = 1L, ...)
{
    derivs <- as.integer(derivs)
    if (derivs < 1L)
        stop("'derivs' has to be a positive integer.")

    degree <- attr(expr, "degree")
    if (derivs > degree + 1L) {
        expr[, ] <- 0
        return(expr)
    }

    if (derivs == 1L)
        return(attr(expr, "bsMat"))

    stop("under developing")
}


##' @rdname deriv
##' @export
deriv.iSpline <- function(expr, derivs = 1L, ...)
{
    derivs <- as.integer(derivs)
    if (derivs < 1L)
        stop("'derivs' has to be a positive integer.")

    degree <- attr(expr, "degree")
    if (derivs > degree + 1L) {
        expr[, ] <- 0
        return(expr)
    }

    if (derivs == 1L)
        return(attr(expr, "msMat"))

    stop("under developing")
}


##' @rdname deriv
##' @export
deriv.cSpline <- function(expr, derivs = 1L, ...)
{
    derivs <- as.integer(derivs)
    if (derivs < 1L)
        stop("'derivs' has to be a positive integer.")

    degree <- attr(expr, "degree")
    if (derivs > degree + 1L) {
        expr[, ] <- 0
        return(expr)
    }

    if (derivs == 1L)
        return(attr(expr, "isMat"))
    if (derivs == 2L)
        return(attr(expr, "msMat"))

    stop("under developing")
}


##' @rdname deriv
deriv.bSpline2 <- function(expr, derivs = 1L, ...) {
    ## FIXME
}


### internal function ==========================================================
derivBs <- function(bsMat, derivs = 1L, ...) {

    ## recover splines information from attributes
    degree <- attr(bsMat, "degree")
    knots <- attr(bsMat, "knots")
    Boundary.knots <- attr(bsMat, "Boundary.knots")
    intercept <- attr(bsMat, "intercept")
    df <- degree + length(knots) + as.integer(intercept)
    x <- attr(bsMat, "x")

    ## for degree zero
    if (derivs > degree) {
        bsMat[, ] <- 0
        return(bsMat)
    }

    ## for degree >= 1L
    ## bSpline()
}
