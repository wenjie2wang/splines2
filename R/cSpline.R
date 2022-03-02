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

##' C-Spline Basis for Polynomial Splines
##'
##' Generates the convex regression spline (called C-spline) basis matrix by
##' integrating I-spline basis for a polynomial spline or the corresponding
##' derivatives.
##'
##' It is an implementation of the closed-form C-spline basis derived from the
##' recursion formula of I-splines and M-splines.
##'
##' @inheritParams bSpline
##'
##' @param degree The degree of C-spline defined to be the degree of the
##'     associated M-spline instead of actual polynomial degree. For example,
##'     C-spline basis of degree 2 is defined as the scaled double integral of
##'     associated M-spline basis of degree 2.
##' @param intercept If \code{TRUE} by default, all of the spline basis
##'     functions are returned.  Notice that when using C-Spline for
##'     shape-restricted regression, \code{intercept = TRUE} should be set even
##'     when an intercept term is considered additional to the spline basis in
##'     the model.
##' @param derivs A nonnegative integer specifying the order of derivatives of
##'     C-splines. The default value is \code{0L} for C-spline basis functions.
##' @param scale A logical value indicating if scaling C-splines is required. If
##'     \code{TRUE} by default, each C-spline basis is scaled to have unit
##'     height at right boundary knot. The corresponding I-spline and M-spline
##'     produced by \code{deriv} methods will be scaled to the same extent.
##'
##' @inherit bSpline return
##'
##' @references
##' Meyer, M. C. (2008). Inference using shape-restricted regression splines.
##' \emph{The Annals of Applied Statistics}, 2(3), 1013--1033.
##'
##' @example inst/examples/ex-cSpline.R
##'
##' @seealso
##' \code{\link{iSpline}} for I-splines;
##' \code{\link{mSpline}} for M-splines.
##'
##' @export
cSpline <- function(x, df = NULL, knots = NULL, degree = 3L,
                    intercept = TRUE, Boundary.knots = NULL,
                    derivs = 0L, scale = TRUE, ...)
{
    ## check inputs
    if ((derivs <- as.integer(derivs)) < 0) {
        stop("The 'derivs' must be a nonnegative integer.")
    }
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
    out <- if (scale) {
               rcpp_cSpline(
                   x = xx,
                   df = df,
                   degree = degree,
                   internal_knots = knots,
                   boundary_knots = Boundary.knots,
                   derivs = derivs,
                   complete_basis = intercept
               )
           } else {
               if (derivs == 0) {
                   rcpp_iSpline(
                       x = xx,
                       df = df,
                       degree = degree,
                       internal_knots = knots,
                       boundary_knots = Boundary.knots,
                       derivs = 0,
                       integral = TRUE,
                       complete_basis = intercept
                   )
               } else {
                   rcpp_iSpline(
                       x = xx,
                       df = df,
                       degree = degree,
                       internal_knots = knots,
                       boundary_knots = Boundary.knots,
                       derivs = derivs - 1,
                       integral = FALSE,
                       complete_basis = intercept
                   )
               }
           }
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
    ## add dimnames for consistency
    name_x <- names(x)
    if (! is.null(name_x)) {
        row.names(out) <- name_x
    }
    ## add class
    if (scale || derivs == 0) {
        ## add "scale" to attributes for predict(), etc.
        attr(out, "scale") <- scale
        class(out) <- c("matrix", "cSpline")
    } else if (derivs == 1) {
        class(out) <- c("matrix", "iSpline")
    } else {
        class(out) <- c("matrix", "mSpline")
    }
    ## return
    out
}
