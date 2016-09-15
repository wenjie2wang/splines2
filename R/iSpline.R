################################################################################
##
##   R package splines2 by Wenjie Wang and Jun Yan
##   Copyright (C) 2016
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


##' I-Spline Basis for Polynomial Splines
##'
##' This function generates the I-spline (integral of M-spline) basis matrix
##' for a polynomial spline.
##'
##' It is an implementation of the close form I-spline basis based on the
##' recursion formula of B-spline basis.  Internally, it calls
##' \code{\link{mSpline}} and \code{\link{bSpline}}, and generates a basis
##' matrix for representing the family of piecewise polynomials and their
##' corresponding integrals with the specified interior knots and degree,
##' evaluated at the values of \code{x}.
##'
##' @param x The predictor variable.  Missing values are allowed and will be
##' returned as they were.
##' @param df Degrees of freedom.  One can specify \code{df} rather than
##' \code{knots}, then the function chooses "df - degree"
##' (minus one if there is an intercept) knots at suitable quantiles of \code{x}
##' (which will ignore missing values).  The default, \code{NULL}, corresponds
##' to no inner knots, i.e., "degree - intercept".
##' @param knots The internal breakpoints that define the spline.  The
##' default is \code{NULL}, which results in a basis for ordinary
##' polynomial regression.  Typical values are the mean or median
##' for one knot, quantiles for more knots.  See also
##' \code{Boundary.knots}.
##' @param degree Non-negative integer degree of the piecewise polynomial. The
##' default value is 3 for cubic splines. Note that the degree of I-spline is
##' defined to be the degree of the associated M-spline instead of actual
##' polynomial degree. In other words, I-spline basis of degree 2 is defined as
##' the integral of associated M-spline basis of degree 2.
##' @param intercept If \code{TRUE}, an intercept is included in the basis;
##' Default is \code{FALSE}.
##' @param Boundary.knots Boundary points at which to anchor the I-spline basis.
##' By default, they are the range of the non-\code{NA} data.  If both
##' \code{knots} and \code{Boundary.knots} are supplied, the basis parameters
##' do not depend on \code{x}. Data can extend beyond \code{Boundary.knots}.
##' @param ... Optional arguments for future usage.
##' @return A matrix of dimension \code{length(x)} by
##' \code{df = degree + length(knots)} (plus on if intercept is included).
##' Attributes that correspond to the arguments specified are returned
##' for usage for \code{\link{predict.iSpline}}. The corresponding M-spline
##' basis matrix is also returned in attribute named \code{msMat}.
##' @references
##' Ramsay, J. O. (1988). Monotone regression splines in action.
##' \emph{Statistical science}, 3(4), 425--441.
##' @examples
##' ## Example given in the reference paper by Ramsay (1988)
##' library(graphics)
##' x <- seq(0, 1, by = .01)
##' knots <- c(0.3, 0.5, 0.6)
##' iMat <- iSpline(x, knots = knots, degree = 2, intercept = TRUE)
##' matplot(x, iMat, type = "l", ylab = "I-spline basis")
##' abline(v = knots, lty = 2, col = "gray")
##' @seealso
##' \code{\link{predict.iSpline}} for evaluation at given (new) values;
##' \code{\link{mSpline}} for M-spline basis;
##' \code{\link{cSpline}} for C-spline basis;
##' \code{\link{bSpline}} for B-spline basis;
##' \code{\link{ibs}} for integral of B-spline basis.
##' @importFrom stats stepfun
##' @export
iSpline <- function(x, df = NULL, knots = NULL, degree = 3, intercept = FALSE,
                    Boundary.knots = range(x), ...) {

    ## M-spline basis for inputs
    msOut <- mSpline(x = x, df = df, knots = knots, degree = degree,
                     intercept = intercept, Boundary.knots = Boundary.knots)

    ## update input
    degree <- attr(msOut, "degree")
    knots <- attr(msOut, "knots")
    bKnots <- attr(msOut, "Boundary.knots")
    ord <- 1L + degree

    ## define knot sequence
    aKnots <- sort(c(rep(bKnots, ord + 1), knots))

    ## function determining j from x
    foo <- stats::stepfun(x = knots, y = seq(ord, length(knots) + ord))
    j <- as.integer(foo(x))

    ## calculate I-spline basis at each x
    ## directly based on B-spline
    bsOut1 <- bSpline(x = x, knots = knots, degree = ord,
                      intercept = FALSE, Boundary.knots = bKnots)
    df <- length(knots) + ord
    bsAugMat <- cbind(j, bsOut1)
    imsOut <- t(apply(bsAugMat, 1, function(b, idx = seq_len(df)) {
        j <- b[1L]
        a <- b[- 1L]
        js <- seq_len(j)
        a[- js] <- 0
        a[js] <- rev(cumsum(rev(a[js])))
        a[idx < j - ord] <- 1            # <=> a[idx < j - degree] <- 1
        a
    }))

    ## Or based on M-spline
    ## generate M-spline basis with (degree + 1)

    ## msOut1 <- mSpline(x = x, knots = knots, degree = ord,
    ##                   intercept = FALSE, Boundary.knots = bKnots)
    ## df <- length(knots) + ord
    ## numer1 <- diff(aKnots, lag = ord + 1)[- 1L]
    ## msMat <- rep(numer1, each = length(x)) * msOut1 / (ord + 1)
    ## msAugMat <- cbind(j, msMat)
    ## imsOut <- t(apply(msAugMat, 1, function(b, idx = seq_len(df)) {
    ##     j <- b[1L]
    ##     a <- b[- 1L]
    ##     js <- seq_len(j)
    ##     a[- js] <- 0
    ##     a[js] <- rev(cumsum(rev(a[js])))
    ##     a[idx < j - ord] <- 1            # <=> a[idx < j - degree] <- 1
    ##     a
    ## }))

    ## intercept
    if (! intercept) imsOut <- imsOut[, - 1L, drop = FALSE]

    ## output
    attributes(imsOut) <- c(attributes(msOut), list(msMat = msOut))
    class(imsOut) <- c("iSpline", "basis", "matrix")
    imsOut
}
