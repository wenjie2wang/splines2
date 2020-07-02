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

##' M-Spline Basis for Polynomial Splines and its Derivatives
##'
##' This function generates the basis matrix of the regression spline called
##' M-spline or its derivatives of given order.  For monotone regression,
##' \code{\link{iSpline}} should be used.
##'
##' It is an implementation of the close form M-spline basis based on
##' relationship between M-spline basis and B-spline basis.  In fact, M-spline
##' basis is a rescaled version of B-spline basis. Internally, it calls function
##' \code{\link{bSpline}} and generates a basis matrix for representing the
##' family of piecewise polynomials with the specified interior knots and
##' degree, evaluated at the values of \code{x}.
##'
##' @inheritParams bSpline
##' @param derivs A non-negative integer specifying the order of derivatives of
##'     M-splines. The default value is \code{0L} for M-spline bases.
##'
##' @return
##' A matrix of dimension \code{length(x)} by
##' \code{df = degree + length(knots)} (plus one if intercept is included).
##' Attributes that correspond to the arguments specified are returned
##' for usage of other functions in this package.
##'
##' @references
##' Ramsay, J. O. (1988). Monotone regression splines in action.
##' \emph{Statistical science}, 3(4), 425--441.
##'
##' @examples
##' ## Example given in the reference paper by Ramsay (1988)
##' library(splines2)
##' x <- seq.int(0, 1, 0.01)
##' knots <- c(0.3, 0.5, 0.6)
##' msMat <- mSpline(x, knots = knots, degree = 2, intercept = TRUE)
##'
##' library(graphics)
##' matplot(x, msMat, type = "l", ylab = "M-spline basis")
##' abline(v = knots, lty = 2, col = "gray")
##'
##' ## derivatives of M-splines
##' dmsMat <- mSpline(x, knots = knots, degree = 2,
##'                   intercept = TRUE, derivs = 1)
##' ## or using the 'deriv' method
##' dmsMat1 <- deriv(msMat)
##' stopifnot(all.equal(dmsMat, dmsMat1, check.attributes = FALSE))
##'
##' @seealso
##' \code{\link{predict.mSpline}} for evaluation at given (new) values;
##' \code{\link{deriv.mSpline}} for derivative method;
##' \code{\link{bSpline}} for B-splines;
##' \code{\link{iSpline}} for I-splines;
##' \code{\link{cSpline}} for C-splines.
##'
##' @export
mSpline <- function(x, df = NULL, knots = NULL, degree = 3L,
                    intercept = FALSE, Boundary.knots = NULL,
                    derivs = 0L, ...)
{
    ## derivs will be checked for range in c++ code
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
