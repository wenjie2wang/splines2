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

##' Derivatives of B-Splines
##'
##' Produces the derivatives of given order of B-splines.
##'
##' This function provides a more user-friendly interface and a more consistent
##' handling for \code{NA}'s than \code{splines::splineDesign()} for derivatives
##' of B-splines.  The implementation is based on the closed-form recursion
##' formula.  At knots, the derivative is defined to be the right derivative
##' except at the right boundary knot.
##'
##' @inheritParams bSpline
##'
##' @param derivs A positive integer specifying the order of derivative.  The
##'     default value is \code{1L} for the first derivative.
##'
##' @inherit bSpline return
##'
##' @references
##' De Boor, Carl. (1978). \emph{A practical guide to splines}.
##' Vol. 27. New York: Springer-Verlag.
##'
##' @example inst/examples/ex-dbs.R
##'
##' @seealso
##' \code{\link{bSpline}} for B-splines;
##' \code{\link{ibs}} for integrals of B-splines.
##'
##' @export
dbs <- function(x, derivs = 1L, df = NULL, knots = NULL, degree = 3L,
                intercept = FALSE, Boundary.knots = NULL, ...)
{
    ## check inputs
    if ((derivs <- as.integer(derivs)) <= 0) {
        stop("The 'derivs' must be a positive integer.")
    }
    if ((degree <- as.integer(degree)) < 0) {
        stop("The 'degree' must be a nonnegative integer.")
    }
    if (is.null(df)) {
        df <- 0L
    } else {
        df <- as.integer(df)
        if (df < 0) {
            stop("The 'df' must be a nonnegative integer.")
        }
    }
    knots <- null2num0(knots)
    Boundary.knots <- null2num0(Boundary.knots)
    ## take care of possible NA's in `x`
    nax <- is.na(x)
    if (all(nax)) {
        stop("The 'x' cannot be all NA's!")
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
            "may cause ill-conditioned basis functions."
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
    class(out) <- c("matrix", "dbs", "splines2")
    ## return
    out
}
