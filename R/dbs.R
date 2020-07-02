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

##' Derivatives of B-Spline Basis
##'
##' Produces the derivative of given order of B-splines.
##'
##' It is an implementation of the close form derivative of B-spline basis based
##' on recursion relation.  At knots, the derivative is defined to be the right
##' derivative.
##'
##' The function is similar with \code{splines::splineDesign}. However, it
##' provides a more user-friendly interface, a more considerate \code{NA}'s
##' handling.
##'
##' @inheritParams bSpline
##' @param derivs A positive integer specifying the order of derivative.  By
##'     default, it is \code{1L} for the first derivative.
##'
##' @return
##' A matrix of dimension \code{length(x)} by
##' \code{df = degree + length(knots)} (plus on if intercept is included).
##' Attributes that correspond to the arguments specified are returned
##' for usage of other functions in this package.
##'
##' @references
##' De Boor, Carl. (1978). \emph{A practical guide to splines}.
##' Vol. 27. New York: Springer-Verlag.
##'
##' @example inst/examples/ex-dbs.R
##'
##' @seealso
##' \code{\link{predict.dbs}} for evaluation at given (new) values;
##' \code{\link{deriv.dbs}} for derivative method;
##' \code{\link{bSpline}} for B-splines;
##' \code{\link{ibs}} for integral of B-splines.
##'
##' @export
dbs <- function(x, derivs = 1L, df = NULL, knots = NULL, degree = 3L,
                intercept = FALSE, Boundary.knots = NULL, ...)
{
    ## derivs will be checked for range in c++ code
    derivs <- as.integer(derivs)
    ## check and reformat 'degree'
    if ((degree <- as.integer(degree)) < 0) {
        stop("'degree' must be a nonnegative integer.")
    }
    if (is.null(df)) {
        df <- 0L
    } else {
        df <- as.integer(df)
        if (df < 0) {
            stop("'df' must be a nonnegative integer.")
        }
    }
    knots <- null2num0(knots)
    Boundary.knots <- null2num0(Boundary.knots)
    ## take care of possible NA's in `x`
    nax <- is.na(x)
    if (all(nax)) {
        stop("'x' cannot be all NA's!")
    }
    nas <- any(nax)
    ## remove NA's
    xx <- if (nas <- any(nax)) {
              x[! nax]
          } else {
              x
          }
    ## call the engine function
    out <- rcpp_bSpline_derivative(
        x = xx,
        derivs = derivs,
        df = df,
        degree = degree,
        internal_knots = knots,
        boundary_knots = Boundary.knots,
        complete_basis = intercept
    )
    ## throw warning if any x is outside of the boundary
    b_knots <- attr(out, "Boundary.knots")
    if (any((xx < b_knots[1L]) | (xx > b_knots[2L]))) {
        warning(wrapMessages(
            "Some 'x' values beyond boundary knots",
            "may cause ill-conditioned bases."
        ))
    }
    ## keep NA's as is
    if (nas) {
        nmat <- matrix(NA, length(nax), ncol(out))
        nmat[! nax, ] <- out
        saved_attr <- attributes(out)
        saved_attr$dim[1] <- length(nax)
        out <- nmat
        attributes(out) <- saved_attr
        attr(out, "x") <- x
    }
    ## add dimnames for consistency with returns from splines::bs
    name_x <- names(x)
    if (! is.null(name_x)) {
        row.names(out) <- name_x
    }
    ## add class
    class(out) <- c("matrix", "dbs")
    ## return
    out
}
