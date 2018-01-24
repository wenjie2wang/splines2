################################################################################
##
##   R package splines2 by Wenjie Wang and Jun Yan
##   Copyright (C) 2016-2018
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


##' C-Spline Basis for Polynomial Splines
##'
##' This function generates the convex regression spline (called C-spline) basis
##' matrix by integrating I-spline basis for a polynomial spline.
##'
##' It is an implementation of the close form C-spline basis derived from
##' the recursion formula of I-spline and M-spline.  Internally, it calls
##' \code{\link{iSpline}} and generates a basis matrix for representing the
##' family of piecewise polynomials and their corresponding integrals with the
##' specified interior knots and degree, evaluated at the values of \code{x}.
##'
##' @usage
##' cSpline(x, df = NULL, knots = NULL, degree = 3L, intercept = FALSE,
##'         Boundary.knots = range(x, na.rm = TRUE), scale = TRUE, ...)
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
##'     default value is 3 for cubic splines.
##' @param intercept If \code{TRUE}, an intercept is included in the basis;
##'     Default is \code{FALSE}.
##' @param Boundary.knots Boundary points at which to anchor the C-spline basis.
##'     By default, they are the range of the non-\code{NA} data.  If both
##'     \code{knots} and \code{Boundary.knots} are supplied, the basis
##'     parameters do not depend on \code{x}. Data can extend beyond
##'     \code{Boundary.knots}.
##' @param scale Logical value (\code{TRUE} by default) indicating whether
##'     scaling on C-spline basis is required. If TRUE, C-spline basis is scaled
##'     to have unit height at right boundary knot; the corresponding I-spline
##'     and M-spline basis matrices shipped in attributes are also scaled to the
##'     same extent.
##' @param ... Optional arguments for future usage.
##'
##' @return A matrix of dimension \code{length(x)} by
##' \code{df = degree + length(knots)} (plus on if intercept is included).
##' The attributes that correspond to the arguments specified are returned
##' for the usage of other functions in this package.
##' @references
##' Meyer, M. C. (2008). Inference using shape-restricted regression splines.
##' \emph{The Annals of Applied Statistics}, 1013--1033. Chicago
##' @examples
##' library(splines2)
##' x <- seq.int(0, 1, 0.01)
##' knots <- c(0.3, 0.5, 0.6)
##'
##' ### when 'scale = TRUE' (by default)
##' csMat <- cSpline(x, knots = knots, degree = 2, intercept = TRUE)
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
##' csMat <- cSpline(x, knots = knots, degree = 2,
##'                  intercept = TRUE, scale = FALSE)
##' ## the corresponding I-splines and M-splines (with same arguments)
##' isMat <- iSpline(x, knots = knots, degree = 2, intercept = TRUE)
##' msMat <- mSpline(x, knots = knots, degree = 2, intercept = TRUE)
##' ## or using deriv methods (much more efficient)
##' isMat1 <- deriv(csMat)
##' msMat1 <- deriv(csMat, derivs = 2)
##' ## equivalent
##' stopifnot(all.equal(isMat, isMat1, check.attributes = FALSE))
##' stopifnot(all.equal(msMat, msMat1, check.attributes = FALSE))
##' @seealso
##' \code{\link{predict.cSpline}} for evaluation at given (new) values;
##' \code{\link{deriv.cSpline}} for derivatives;
##' \code{\link{iSpline}} for I-splines;
##' \code{\link{mSpline}} for M-splines.
##' @importFrom stats stepfun
##' @export
cSpline <- function(x, df = NULL, knots = NULL, degree = 3L, intercept = FALSE,
                    Boundary.knots = range(x, na.rm = TRUE), scale = TRUE, ...)
{
    ## I-spline basis for inputs
    isOut <- iSpline(x = x, df = df, knots = knots, degree = degree,
                     intercept = intercept, Boundary.knots = Boundary.knots)

    ## update input
    degree <- attr(isOut, "degree")
    knots <- attr(isOut, "knots")
    bKnots <- attr(isOut, "Boundary.knots")
    ord <- 1L + degree
    nKnots <- length(knots)
    df <- nKnots + ord

    ## take care of possible NA's in `x` for the following calculation
    nax <- is.na(x)
    if ((nas <- any(nax)))
        x <- x[! nax]
    nX <- length(x)

    ## define knot sequence
    aKnots <- sort(c(rep(bKnots, ord + 1L), knots))

    ## generate I-spline basis with (degree + 1)
    augX <- c(x, bKnots[2L])
    isOut1 <- iSpline(x = augX, knots = knots, degree = ord,
                      intercept = FALSE, Boundary.knots = bKnots)

    ## function determining j from x
    j <- if (length(knots)) {
             foo <- stats::stepfun(x = knots, y = seq.int(ord, df))
             as.integer(foo(augX))
         } else {
             rep.int(ord, nX + 1L)
         }

    numer1 <- diff(aKnots, lag = ord + 1L)[- 1L]
    ## if there is at least one internal knot
    if (nKnots) {
        ## calculate C-spline basis at each internal knot t_j
        isOutKnots <- iSpline(knots, knots = knots, degree = ord,
                              intercept = FALSE, Boundary.knots = bKnots)
        matKnots <- rep(numer1, each = nKnots) * isOutKnots / (ord + 1)
        augKnots <- seq_len(nKnots) + ord
        diffKnots <- diff(knots)
        csKnots <- lapply(seq_len(nKnots), function(i, idx) {
            ji <- augKnots[i]
            a <- matKnots[i, ]
            js <- seq_len(ji)
            a[- js] <- 0
            a[js] <- rev(cumsum(rev(a[js])))
            a[idx < ji - ord] <- diffKnots[ji - ord - 1L]
            a
        }, idx = seq_len(df))
        csKnots <- do.call(rbind, csKnots)

        idxMat <- lower.tri(csKnots, diag = TRUE)
        linList <- lapply(seq_len(nKnots), function(ind) {
            cumsum(csKnots[idxMat[, ind], ind])
        })
        csKnots[idxMat] <- do.call(c, linList)
    } else {
        csKnots <- matrix(0, 1L, df)
    }

    ## calculate C-spline basis at each x
    matX <- rep(numer1, each = nX + 1) * isOut1 / (ord + 1)
    csOut <- lapply(seq_len(nX + 1L), function(i, idx) {
        ji <- j[i]
        xx <- augX[i]
        a <- matX[i, ]
        js <- seq_len(ji)
        a[- js] <- 0
        a[js] <- rev(cumsum(rev(a[js])))
        a[idx < ji - ord] <- xx - knots[ji - ord] +
            csKnots[ji - ord, idx < ji - ord]
        a
    }, idx = seq_len(df))
    csOut <- do.call(rbind, csOut)

    if (! intercept)
        csOut <- csOut[, - 1L, drop = FALSE]
    scl <- unname(csOut[nX + 1L, ])
    csOut <- csOut[- (nX + 1L), ]

    ## mSpline basis matrix
    msMat <- attr(isOut, "msMat")

    ## keep NA's as is for csOut
    if (nas) {
        nmat <- matrix(NA, length(nax), ncol(csOut))
        nmat[! nax, ] <- csOut
        csOut <- nmat
    }

    ## scale C-spline, I-spline, and M-spline basis
    if (scale) {
        vec <- rep(1 / scl, each = length(nax))
        csOut <- vec * csOut
        isOut <- vec * isOut
        msMat <- vec * msMat
        attr(isOut, "scale") <- attr(msMat, "scale") <- scale
        attr(isOut, "scales") <- attr(msMat, "scales") <- scl
    }

    ## output
    attr(isOut, "msMat") <- NULL
    attributes(csOut) <- c(attributes(isOut),
                           list(isMat = isOut, msMat = msMat,
                                scale = scale, scales = scl))
    attr(csOut, "derivs") <- NULL
    class(csOut) <- c("matrix", "cSpline")
    csOut
}
