################################################################################
##
##   R package splines2 by Wenjie Wang and Jun Yan
##   Copyright (C) 2016-2019
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
##' @usage
##' mSpline(x, df = NULL, knots = NULL, degree = 3L, intercept = FALSE,
##'         Boundary.knots = range(x, na.rm = TRUE), derivs = 0L, ...)
##'
##' @param x The predictor variable.  Missing values are allowed and will be
##'     returned as they were.
##' @param df Degrees of freedom.  One can specify \code{df} rather than
##'     \code{knots}, then the function chooses "df - degree" (minus one if
##'     there is an intercept) knots at suitable quantiles of \code{x} (which
##'     will ignore missing values).  The default, \code{NULL}, corresponds to
##'     no inner knots, i.e., "degree - intercept".
##' @param knots The internal breakpoints that define the spline.  The default
##'     is \code{NULL}, which results in a basis for ordinary polynomial
##'     regression.  Typical values are the mean or median for one knot,
##'     quantiles for more knots.  See also \code{Boundary.knots}.
##' @param degree Non-negative integer degree of the piecewise polynomial. The
##'     default value is 3 for cubic splines. Zero degree is allowed for
##'     piecewise constant basis.
##' @param intercept If \code{TRUE}, all bases will be returned. The default
##'     value is \code{FALSE}.
##' @param Boundary.knots Boundary points at which to anchor the M-spline basis.
##'     By default, they are the range of the non-\code{NA} data.  If both
##'     \code{knots} and \code{Boundary.knots} are supplied, the basis
##'     parameters do not depend on \code{x}. Data can extend beyond
##'     \code{Boundary.knots}.
##' @param derivs A non-negative integer specifying the order of derivatives of
##'     M-splines. The default value is \code{0L} for M-spline bases.
##' @param ... Optional arguments for future usage.
##'
##' @return A matrix of dimension \code{length(x)} by
##' \code{df = degree + length(knots)} (plus one if intercept is included).
##' Attributes that correspond to the arguments specified are returned
##' for usage of other functions in this package.
##' @references
##' Ramsay, J. O. (1988). Monotone regression splines in action.
##' \emph{Statistical science}, 3(4), 425--441.
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
##' @seealso
##' \code{\link{predict.mSpline}} for evaluation at given (new) values;
##' \code{\link{deriv.mSpline}} for derivative method;
##' \code{\link{bSpline}} for B-splines;
##' \code{\link{iSpline}} for I-splines;
##' \code{\link{cSpline}} for C-splines.
##' @export
mSpline <- function(x, df = NULL, knots = NULL, degree = 3L, intercept = FALSE,
                    Boundary.knots = range(x, na.rm = TRUE), derivs = 0L, ...)
{
    ## check order of derivative
    if (! missing(derivs)) {
        derivs <- as.integer(derivs)
        if (derivs < 0L)
            stop("'derivs' has to be a non-negative integer.")
    }

    bsOut <- if (derivs) {
                 dbs(x = x, derivs = derivs, df = df, knots = knots,
                     degree = degree, intercept = intercept,
                     Boundary.knots = Boundary.knots, ...)
             } else {
                 bSpline(x = x, df = df, knots = knots,
                         degree = degree, intercept = intercept,
                         Boundary.knots = Boundary.knots, ...)
             }

    ## update input
    ord <- attr(bsOut, "degree") + 1L
    knots <- attr(bsOut, "knots")
    bKnots <- attr(bsOut, "Boundary.knots")

    ## define knot sequence
    aKnots <- sort(c(rep(bKnots, ord), knots))

    ## transformation from B-splines to M-splines
    denom <- diff(aKnots, lag = ord)
    transCoef <- ifelse(denom > 0, ord / denom, 0)
    if (! intercept)
        transCoef <- transCoef[- 1L]
    msOut <- rep(transCoef, each = length(x)) * bsOut
    attr(msOut, "derivs") <- derivs
    class(msOut) <- c("matrix", "mSpline")
    msOut
}
