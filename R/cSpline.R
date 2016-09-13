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


##' C-Spline Basis for Polynomial Splines
##'
##' This function generates the C-spline (integral of I-spline) basis matrix
##' for a polynomial spline.
##'
##' It is an implementation of the close form C-spline basis based
##' on recursion relation.  Internally, it calls \code{\link{iSpline}} and
##' generates a basis matrix for representing the family of piecewise
##' polynomials and their corresponding integrals with the specified interior
##' knots and degree, evaluated at the values of \code{x}.
##' @param x The predictor variable.  Missing values are allowed and will be
##' returned but ignored in computation.
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
##' @param degree Degree of the piecewise polynomial. The default value is 3
##' for cubic splines.
##' @param intercept If \code{TRUE}, an intercept is included in the basis;
##' Default is \code{FALSE}.
##' @param Boundary.knots Boundary points at which to anchor the C-spline basis.
##' By default, they are the range of the non-\code{NA} data.  If both
##' \code{knots} and \code{Boundary.knots} are supplied, the basis parameters
##' do not depend on \code{x}. Data can extend beyond \code{Boundary.knots}.
##' @param rescale Logical value (\code{TRUE} by default) indicating whether the
##' rescaling on C-spline basis is required. If TRUE, C-spline basis is rescaled
##' to have unit height at right boundary knot; the corresponding I-spline and
##' M-spline basis matrices are also rescaled to the same extent.
##' @param ... Optional arguments for future usage.
##' @return A matrix of dimension \code{length(x)} by
##' \code{df = degree + length(knots)} (plus on if intercept is included).
##' Attributes that correspond to the arguments specified are returned
##' for usage for \code{\link{predict.iSpline}}. The corresponding M-spline
##' basis matrix is also returned in attribute named \code{msMat}.
##' @references
##' FIXME
##' @examples
##' x <- seq(0, 1, by = .01)
##' knots <- c(0.3, 0.5, 0.6)
##' cMat <- cSpline(x, knots = knots, degree = 2, intercept = TRUE)
##' matplot(x, cMat, type = "l", ylab = "C-spline basis")
##' abline(v = knots, lty = 2, col = "gray")
##' @seealso
##' \code{\link{predict.cSpline}} for evaluation at given (new) values;
##' \code{\link{iSpline}} for I-spine basis;
##' \code{\link{mSpline}} for M-spine basis;
##' \code{\link{ibs}} for integral of B-spline basis.
##' @importFrom stats stepfun
##' @export
cSpline <- function(x, df = NULL, knots = NULL, degree = 3, intercept = FALSE,
                   Boundary.knots = range(x), rescale = TRUE, ...) {

    ## I-spline basis for inputs
    imsOut <- iSpline(x = x, df = df, knots = knots, degree = degree,
                     intercept = intercept, Boundary.knots = Boundary.knots)

    ## update input
    degree <- attr(imsOut, "degree")
    knots <- attr(imsOut, "knots")
    bKnots <- attr(imsOut, "Boundary.knots")
    ord <- 1L + degree
    nKnots <- length(knots)
    df <- nKnots + ord
    nX <- length(x)

    ## define knot sequence
    aKnots <- sort(c(rep(bKnots, ord + 1), knots))

    ## generate I-spline basis with (degree + 1)
    augX <- c(x, bKnots[2])
    imsOut1 <- iSpline(x = augX, knots = knots, degree = ord,
                      intercept = FALSE, Boundary.knots = bKnots)

    ## function determining j from x
    foo <- stats::stepfun(x = knots, y = seq(ord, df))
    j <- as.integer(foo(augX))

    ## calculate C-spline basis at each internal knot t_j
    numer1 <- diff(aKnots, lag = ord + 1)[- 1L]
    imsOutKnots <- iSpline(knots, knots = knots, degree = ord,
                          intercept = FALSE, Boundary.knots = bKnots)
    matKnots <- rep(numer1, each = nKnots) * imsOutKnots / (ord + 1)
    augMatKnots <- cbind(seq_len(nKnots) + ord, matKnots)
    diffKnots <- diff(knots)
    csKnots <- t(apply(augMatKnots, 1, function(b, idx = seq_len(df)) {
        j <- b[1L]
        a <- b[- 1L]
        js <- seq_len(j)
        a[- js] <- 0
        a[js] <- rev(cumsum(rev(a[js])))
        a[idx < j - ord] <- diffKnots[j - ord - 1]
        a
    }))
    idxMat <- lower.tri(csKnots, diag = TRUE)
    linList <- lapply(seq_len(nKnots), function (ind) {
        cumsum(csKnots[idxMat[, ind], ind])
    })
    csKnots[idxMat] <- do.call("c", linList)

    ## calculate C-spline basis at each x
    matX <- rep(numer1, each = nX + 1) * imsOut1 / (ord + 1)
    augMatX <- cbind(j, augX, matX)
    csOut <- t(apply(augMatX, 1, function(b, idx = seq_len(df)) {
        j <- b[1L]
        xx <- b[2L]
        a <- b[- seq_len(2)]
        js <- seq_len(j)
        a[- js] <- 0
        a[js] <- rev(cumsum(rev(a[js])))
        a[idx < j - ord] <- xx - knots[j - ord] +
            csKnots[j - ord, idx < j - ord]
        a
    }))
    if (! intercept) csOut <- csOut[, - 1L, drop = FALSE]
    scl <- csOut[nX + 1, ]
    csOut <- csOut[- (nX + 1), ]

    ## mSpline basis matrix
    msMat <- attr(imsOut, "msMat")
    ## rescale C-spline, I-spline, and M-spline basis
    if (rescale) {
        vec <- rep(1 / scl, each = nX)
        csOut <- vec * csOut
        imsOut <- vec * imsOut
        msMat <- vec * msMat
    }

    ## output
    attr(imsOut, "msMat") <- NULL
    attributes(csOut) <- c(attributes(imsOut),
                          list(imsMat = imsOut, msMat = msMat,
                               rescale = rescale))
    class(csOut) <- c("cSpline", "basis", "matrix")
    csOut
}
