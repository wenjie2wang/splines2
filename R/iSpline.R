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

##' I-Spline Basis for Polynomial Splines
##'
##' Generates the I-spline (integral of M-spline) basis matrix for a polynomial
##' spline or the corresponding derivatives of given order.
##'
##' It is an implementation of the close form I-spline basis based on the
##' recursion formula given by Ramsay (1988).
##'
##' @inheritParams bSpline
##'
##' @param degree The degree of I-spline defined to be the degree of the
##'     associated M-spline instead of actual polynomial degree. For example,
##'     I-spline basis of degree 2 is defined as the integral of associated
##'     M-spline basis of degree 2.
##' @param intercept If \code{TRUE} by default, all spline bases are included.
##'     Notice that when using I-Spline for monotonic regression,
##'     \code{intercept = TRUE} should be set even when an intercept term is
##'     considered additional to the spline bases in the model.
##' @param derivs A non-negative integer specifying the order of derivatives of
##'     I-splines.
##'
##' @inherit bSpline return
##'
##' @references
##' Ramsay, J. O. (1988). Monotone regression splines in action.
##' \emph{Statistical science}, 3(4), 425--441.
##'
##' @example inst/examples/ex-iSpline.R
##'
##' @seealso
##' \code{\link{mSpline}} for M-splines;
##' \code{\link{cSpline}} for C-splines;
##'
##' @export
iSpline <- function(x, df = NULL, knots = NULL, degree = 3L,
                    intercept = TRUE, Boundary.knots = NULL,
                    derivs = 0L, ...)
{
    ## check inputs
    if ((derivs <- as.integer(derivs)) < 0) {
        stop("The 'derivs' must be a non-negative integer.")
    }
    if (derivs > 0) {
        return(mSpline(x = x,
                       df = df,
                       knots = knots,
                       degree = degree,
                       intercept = intercept,
                       Boundary.knots = Boundary.knots,
                       derivs = derivs - 1L, ...))
    }
    ## else I-Spline basis
    if ((degree <- as.integer(degree)) < 0)
        stop("The 'degree' must be a nonnegative integer.")
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
    ## remove NA's in x
    xx <- if (nas <- any(nax)) {
              x[! nax]
          } else {
              x
          }
    ## call the engine function
    out <- rcpp_iSpline(
        x = xx,
        df = df,
        degree = degree,
        internal_knots = knots,
        boundary_knots = Boundary.knots,
        derivs = derivs,
        integral = FALSE,
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
    ## add dimnames for consistency
    name_x <- names(x)
    if (! is.null(name_x)) {
        row.names(out) <- name_x
    }
    ## add class
    class(out) <- c("matrix", "iSpline")
    out
}
