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

##' M-Spline Basis for Polynomial Splines
##'
##' Generates the basis matrix of the regression spline called M-spline or the
##' corresponding derivatives of given order.  For monotone regression,
##' \code{\link{iSpline}} should be used instead of M-splines.
##'
##' It is an implementation of the close form M-spline basis based on the
##' recursion formula given by Ramsay (1988).
##'
##' @inheritParams bSpline
##'
##' @param derivs A non-negative integer specifying the order of derivatives of
##'     M-splines. The default value is \code{0L} for M-spline bases.
##'
##' @inherit bSpline return
##'
##' @references
##' Ramsay, J. O. (1988). Monotone regression splines in action.
##' \emph{Statistical science}, 3(4), 425--441.
##'
##' @example inst/examples/ex-mSpline.R
##'
##' @seealso
##' \code{\link{bSpline}} for B-splines;
##' \code{\link{iSpline}} for I-splines;
##' \code{\link{cSpline}} for C-splines.
##'
##' @export
mSpline <- function(x, df = NULL, knots = NULL, degree = 3L,
                    intercept = FALSE, Boundary.knots = NULL,
                    derivs = 0L, ...)
{
    ## check inputs
    if ((derivs <- as.integer(derivs)) < 0) {
        stop("'derivs' must be a non-negative integer.")
    }
    if ((degree <- as.integer(degree)) < 0)
        stop("'degree' must be a nonnegative integer.")
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
        stop("The 'x' cannot be all NA's!")
    }
    ## remove NA's in x
    xx <- if (nas <- any(nax)) {
              x[! nax]
          } else {
              x
          }
    out <- if (derivs > 0) {
               ## call the engine function
               rcpp_mSpline_derivative(
                   x = xx,
                   derivs = derivs,
                   df = df,
                   degree = degree,
                   internal_knots = knots,
                   boundary_knots = Boundary.knots,
                   complete_basis = intercept
               )
           } else {
               ## call the engine function
               rcpp_mSpline_basis(
                   x = xx,
                   df = df,
                   degree = degree,
                   internal_knots = knots,
                   boundary_knots = Boundary.knots,
                   complete_basis = intercept
               )
           }
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
    class(out) <- c("matrix", "mSpline")
    out
}
