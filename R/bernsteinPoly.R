##
## R package splines2 by Wenjie Wang and Jun Yan
## Copyright (C) 2016-2022
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

##' Generalized Bernstein Polynomial Basis
##'
##' Returns a generalized Bernstein polynomial basis matrix of given degree over
##' a specified range.
##'
##' @name bernsteinPoly
##'
##' @inheritParams bSpline
##'
##' @param x The predictor variable taking values inside of the specified
##'     boundary.  Missing values are allowed and will be returned as they are.
##' @param integral A logical value.  If \code{TRUE}, the integrals of the
##'     Bernstein polynomials will be returned.  The default value is
##'     \code{FALSE}.
##' @param Boundary.knots Boundary points at which to anchor the Bernstein
##'     polynomial basis. The default value is \code{NULL} and the boundary
##'     knots is set internally to be \code{range(x, na.rm = TRUE)}.
##' @param degree A nonnegative integer representing the degree of the
##'     polynomials.
##' @param derivs A nonnegative integer specifying the order of derivatives.
##'     The default value is \code{0L} for Bernstein polynomial basis functions.
##'
##' @return A numeric matrix of dimension \code{length(x)} by \code{degree +
##'     as.integer(intercept)}.
##'
##' @example inst/examples/ex-bernsteinPoly.R
##'
##' @export
bernsteinPoly <- function(x, degree = 3, intercept = FALSE,
                          Boundary.knots = NULL,
                          derivs = 0L, integral = FALSE,
                          ...)
{
    ## check inputs
    if ((degree <- as.integer(degree)) < 0)
        stop("The 'degree' must be a nonnegative integer.")
    ## take care of possible NA's in `x`
    nax <- is.na(x)
    if (all(nax)) {
        stop("The 'x' cannot be all NA's!")
    }
    ## remove NA's in x
    xx <- if (nas <- any(nax)) {
              x[! nax]
          } else {
              x
          }
    Boundary.knots <- null2num0(Boundary.knots)
    ## call the engine function
    out <- rcpp_bernsteinPoly(
        x = xx,
        degree = degree,
        derivs = derivs,
        integral = integral,
        boundary_knots = Boundary.knots,
        complete_basis = intercept
    )
    ## keep NA's as is
    if (nas) {
        nmat <- matrix(NA_real_, length(nax), ncol(out))
        nmat[! nax, ] <- out
        saved_attr <- attributes(out)
        saved_attr$dim[1] <- length(nax)
        out <- nmat
        attributes(out) <- saved_attr
        attr(out, "x") <- x
    }
    ## add dimnames for consistency with bs returns
    name_x <- names(x)
    if (! is.null(name_x)) {
        row.names(out) <- name_x
    }
    ## add class
    class(out) <- c("matrix", "bernsteinPoly")
    ## return
    out
}
