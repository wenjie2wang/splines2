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

##' C-Spline Basis for Polynomial Splines or its derivatives
##'
##' Generates the convex regression spline (called C-spline) basis
##' matrix by integrating I-spline basis for a polynomial spline.
##'
##' It is an implementation of the close form C-spline basis derived from the
##' recursion formula of I-spline and M-spline.
##'
##' @inheritParams bSpline
##'
##' @param degree Non-negative integer degree of the piecewise polynomial. The
##'     default value is 3 for cubic splines.
##' @param intercept If \code{TRUE} by default, all spline bases are included.
##'     Notice that when using C-Spline for shape-restricted regression,
##'     \code{intercept = TRUE} should be set even when an intercept term is
##'     considered additional to the spline bases in the model.
##' @param derivs A non-negative integer specifying the order of derivatives of
##'     C-splines. The default value is \code{0L} for C-spline bases.
##' @param scale Logical value (\code{TRUE} by default) indicating whether
##'     scaling on C-spline basis is required. If TRUE, C-spline basis is scaled
##'     to have unit height at right boundary knot; the corresponding I-spline
##'     and M-spline basis matrices shipped in attributes are also scaled to the
##'     same extent.
##'
##' @return
##' A matrix of dimension \code{length(x)} by
##' \code{df = degree + length(knots)} (plus on if intercept is included).
##' The attributes that correspond to the arguments specified are returned
##' for the usage of other functions in this package.
##'
##' @references
##' Meyer, M. C. (2008). Inference using shape-restricted regression splines.
##' \emph{The Annals of Applied Statistics}, 1013--1033. Chicago
##'
##' @examples
##' library(splines2)
##' x <- seq.int(0, 1, 0.01)
##' knots <- c(0.3, 0.5, 0.6)
##'
##' ### when 'scale = TRUE' (by default)
##' csMat <- cSpline(x, knots = knots, degree = 2)
##'
##' library(graphics)
##' matplot(x, csMat, type = "l", ylab = "C-spline basis")
##' abline(v = knots, lty = 2, col = "gray")
##' isMat <- deriv(csMat)
##' msMat <- deriv(csMat, derivs = 2)
##' matplot(x, isMat, type = "l", ylab = "scaled I-spline basis")
##' matplot(x, msMat, type = "l", ylab = "scaled M-spline basis")
##'
##' ### when 'scale = FALSE'
##' csMat <- cSpline(x, knots = knots, degree = 2, scale = FALSE)
##' ## the corresponding I-splines and M-splines (with same arguments)
##' isMat <- iSpline(x, knots = knots, degree = 2)
##' msMat <- mSpline(x, knots = knots, degree = 2, intercept = TRUE)
##' ## or using deriv methods (more efficient)
##' isMat1 <- deriv(csMat)
##' msMat1 <- deriv(csMat, derivs = 2)
##' ## equivalent
##' stopifnot(all.equal(isMat, isMat1, check.attributes = FALSE))
##' stopifnot(all.equal(msMat, msMat1, check.attributes = FALSE))
##'
##' @seealso
##' \code{\link{predict.cSpline}} for evaluation at given (new) values;
##' \code{\link{deriv.cSpline}} for derivatives;
##' \code{\link{iSpline}} for I-splines;
##' \code{\link{mSpline}} for M-splines.
##'
##' @export
cSpline <- function(x, df = NULL, knots = NULL, degree = 3L,
                    intercept = TRUE, Boundary.knots = NULL,
                    derivs = 0L, scale = TRUE, ...)
{
    derivs <- as.integer(derivs)
    ## check inputs
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
    out <- if (scale) {
               if (derivs > 0) {
                   rcpp_cSpline_derivative(
                       x = xx,
                       derivs = derivs,
                       df = df,
                       degree = degree,
                       internal_knots = knots,
                       boundary_knots = Boundary.knots,
                       complete_basis = intercept
                   )
               } else {
                   rcpp_cSpline_basis(
                       x = xx,
                       df = df,
                       degree = degree,
                       internal_knots = knots,
                       boundary_knots = Boundary.knots,
                       complete_basis = intercept
                   )
               }
           } else {
               if (derivs == 0) {
                   rcpp_iSpline_integral(
                       x = xx,
                       df = df,
                       degree = degree,
                       internal_knots = knots,
                       boundary_knots = Boundary.knots,
                       complete_basis = intercept
                   )
               } else if (derivs == 1) {
                   rcpp_iSpline_basis(
                       x = xx,
                       df = df,
                       degree = degree,
                       internal_knots = knots,
                       boundary_knots = Boundary.knots,
                       complete_basis = intercept
                   )
               } else if (derivs > 1) {
                   rcpp_iSpline_derivative(
                       x = xx,
                       derivs = derivs - 1L,
                       df = df,
                       degree = degree,
                       internal_knots = knots,
                       boundary_knots = Boundary.knots,
                       complete_basis = intercept
                   )
               }
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
    if (scale || derivs == 0) {
        class(out) <- c("matrix", "cSpline")
    } else if (derivs == 1) {
        class(out) <- c("matrix", "iSpline")
    } else {
        class(out) <- c("matrix", "mSpline")
    }
    ## return
    out
}
