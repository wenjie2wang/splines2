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


##' Integral of B-Spline Basis for Polynomial Splines
##'
##' This function generates the integral of B-spline basis matrix
##' for a polynomial spline. The arguments are exactly the same with function
##' \code{\link[splines]{bs}} in package \code{splines}.
##'
##' It is an implementation of the close form integral of B-spline basis based
##' on recursion relation.  Internally, it calls \code{\link{bSpline}} and
##' generates a basis matrix for representing the family of piecewise
##' polynomials and their corresponding integrals with the specified interior
##' knots and degree, evaluated at the values of \code{x}.
##'
##' @usage
##' ibs(x, df = NULL, knots = NULL, degree = 3, intercept = FALSE,
##'     Boundary.knots = range(x, na.rm = TRUE), ...)
##'
##' @param x The predictor variable.  Missing values are allowed and will be
##' returned as they were.
##' @param df Degrees of freedom of the B-spline basis to be integrated.
##' One can specify \code{df} rather than \code{knots}, then the function
##' chooses "df - degree" (minus one if there is an intercept) knots at
##' suitable quantiles of \code{x} (which will ignore missing values).
##' The default, \code{NULL}, corresponds to no inner knots, i.e.,
##' "degree - intercept".
##' @param knots The internal breakpoints that define the B-spline basis to be
##' integrated.  The default is \code{NULL}, which results in a basis for
##' ordinary polynomial regression.  Typical values are the mean or median
##' for one knot, quantiles for more knots.  See also \code{Boundary.knots}.
##' @param degree Non-negative integer degree of the piecewise polynomial to be
##' integrated. The default value is 3 for the integral of cubic B-splines.
##' @param intercept If \code{TRUE}, an intercept is included in the basis;
##' Default is \code{FALSE}.
##' @param Boundary.knots Boundary points at which to anchor the B-spline basis
##' to be integrated. By default, they are the range of the non-\code{NA} data.
##' If both \code{knots} and \code{Boundary.knots} are supplied, the basis
##' parameters do not depend on \code{x}. Data can extend beyond
##' \code{Boundary.knots}.
##' @param ... Optional arguments for future usage.
##'
##' @return A matrix of dimension \code{length(x)} by
##' \code{df = degree + length(knots)} (plus on if intercept is included).
##' Attributes that correspond to the arguments specified are returned
##' for usage of other functions in this package.
##' @references
##' De Boor, Carl. (1978). \emph{A practical guide to splines}.
##' Vol. 27. New York: Springer-Verlag.
##' @examples
##' library(splines2)
##' x <- seq(0, 1, 0.01)
##' knots <- c(0.2, 0.4, 0.7, 0.9)
##' ibsMat <- ibs(x, knots = knots, degree = 1, intercept = TRUE)
##'
##' ## obtain the B-spline bases integrated by function bSpline
##' bsMat0 <- bSpline(x, knots = knots, degree = 1, intercept = TRUE)
##' ## or by function deriv (recommended) that directly extracts the existing
##' ## result from the attribute of ibsMat and thus is more efficient.
##' bsMat <- deriv(ibsMat)
##' all.equal(bsMat0, bsMat)                # equivalent
##'
##' ## plot B-spline basis with their corresponding integrals
##' library(graphics)
##' par(mfrow = c(1, 2))
##' matplot(x, bsMat, type = "l", ylab = "B-spline basis")
##' abline(v = knots, lty = 2, col = "gray")
##' matplot(x, ibsMat, type = "l", ylab = "Integral of B-spline basis")
##' abline(v = knots, lty = 2, col = "gray")
##' par(mfrow = c(1, 1))
##' @seealso
##' \code{\link{bSpline}} for B-spline basis;
##' \code{\link{mSpline}} for M-spline basis;
##' \code{\link{iSpline}} for I-spline basis.
##' \code{\link{cSpline}} for C-spline basis.
##' @export
ibs <- function(x, df = NULL, knots = NULL, degree = 3, intercept = FALSE,
                Boundary.knots = range(x, na.rm = TRUE), ...)
{
    ## B-spline basis for inputs
    bsOut <- bSpline(x = x, df = df, knots = knots, degree = degree,
                     intercept = intercept, Boundary.knots = Boundary.knots)

    ## update input
    degree <- attr(bsOut, "degree")
    knots <- attr(bsOut, "knots")
    bKnots <- attr(bsOut, "Boundary.knots")
    ord <- 1L + degree

    ## define knot sequence
    aKnots <- sort(c(rep(bKnots, ord), knots))

    ## generate B-spline basis with (degree + 1)
    bsOut1 <- bSpline(x = x, knots = knots, degree = ord,
                      intercept = FALSE, Boundary.knots = bKnots)
    numer1 <- diff(aKnots, lag = ord)
    if (! intercept) {
        bsOut1 <- bsOut1[, - 1L, drop = FALSE]
        numer1 <- numer1[- 1L]
    }
    numer2 <- apply(bsOut1, 1, function(a) rev(cumsum(rev(a))))
    ibsOut <- t(numer1 * numer2) / ord

    ## output
    attributes(ibsOut) <- c(attributes(bsOut),
                            list(bsMat = bsOut, x = x))
    class(ibsOut) <- c("matrix", "ibs")
    ibsOut
}
