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

##' I-Spline Basis for Polynomial Splines or its derivatives
##'
##' Generates the I-spline (integral of M-spline) basis matrix for a polynomial
##' spline or its derivatives of given order.
##'
##' It is an implementation of the close form I-spline basis based on the
##' recursion formula of B-spline basis.  Internally, it calls
##' \code{\link{mSpline}} and \code{\link{bSpline}}, and generates a basis
##' matrix for representing the family of piecewise polynomials and their
##' corresponding integrals with the specified interior knots and degree,
##' evaluated at the values of \code{x}.
##'
##' @inheritParams bSpline
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
##' @return
##' A matrix of dimension \code{length(x)} by
##' \code{df = degree + length(knots)} (plus on if intercept is included).
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
##' x <- seq.int(0, 1, by = 0.01)
##' knots <- c(0.3, 0.5, 0.6)
##' isMat <- iSpline(x, knots = knots, degree = 2)
##'
##' library(graphics)
##' matplot(x, isMat, type = "l", ylab = "I-spline basis")
##' abline(v = knots, lty = 2, col = "gray")
##'
##' ## the derivative of I-splines is M-spline
##' msMat1 <- iSpline(x, knots = knots, degree = 2, derivs = 1)
##' msMat2 <- mSpline(x, knots = knots, degree = 2, intercept = TRUE)
##' stopifnot(all.equal(msMat1, msMat2))
##'
##' @seealso
##' \code{\link{predict.iSpline}} for evaluation at given (new) values;
##' \code{\link{deriv.iSpline}} for derivative method;
##' \code{\link{mSpline}} for M-splines;
##' \code{\link{cSpline}} for C-splines;
##'
##' @export
iSpline <- function(x, df = NULL, knots = NULL, degree = 3L,
                    intercept = TRUE, Boundary.knots = NULL,
                    derivs = 0L, ...)
{
    derivs <- as.integer(derivs)
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
    ## call the engine function
    out <- rcpp_iSpline_basis(
        x = xx,
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
    ## add dimnames for consistency
    name_x <- names(x)
    if (! is.null(name_x)) {
        row.names(out) <- name_x
    }
    ## add class
    class(out) <- c("matrix", "iSpline")
    out
}
