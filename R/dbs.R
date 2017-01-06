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


##' Derivative of B-Spline Basis for Polynomial Splines
##'
##' This function produces the derivative of B-splines.
##'
##' It is an implementation of the close form derivative of B-spline basis
##' based on recursion relation.  Internally, it calls \code{\link{bSpline}}
##' and generates a basis matrix for representing the family of piecewise
##' polynomials and their corresponding derivative with the specified
##' interior knots and degree, evaluated at the values of \code{x}. The
##' function \code{splineDesign} in \code{splines} package can also be used to
##' calculate derivative of B-splines.
##'
##' @usage
##' dbs(x, derivs = 1L, df = NULL, knots = NULL, degree = 3L,
##'     intercept = FALSE, Boundary.knots = range(x, na.rm = TRUE), ...)
##'
##' @param x The predictor variable.  Missing values are allowed and will be
##' returned as they were.
##' @param derivs A positive integer specifying the order of derivative.
##' By default, it is \code{1L} for the first derivative.
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
##' knots <- c(0.2, 0.4, 0.7)
##' ## the second derivative of cubic B-splines with three internal knots
##' dMat <- dbs(x, derivs = 2L, knots = knots, intercept = TRUE)
##'
##' ## compare with the results from package splines
##' ord <- attr(dMat, "degree") + 1L
##' bKnots <- attr(dMat, "Boundary.knots")
##' aKnots <- c(rep(bKnots[1L], ord), knots, rep(bKnots[2L], ord))
##' res <- splines::splineDesign(aKnots, x = x, derivs = 2L)
##' all.equal(res, dMat, check.attributes = FALSE)
##' @seealso
##' \code{\link{bSpline}} for B-spline basis;
##' \code{\link{mSpline}} for M-spline basis;
##' \code{\link{iSpline}} for I-spline basis.
##' \code{\link{cSpline}} for C-spline basis.
##' @export
dbs <- function(x, derivs = 1L, df = NULL, knots = NULL, degree = 3L,
                intercept = FALSE, Boundary.knots = range(x, na.rm = TRUE), ...)
{
    ## test order of derivative
    derivs <- as.integer(derivs)
    if (derivs < 1L)
        stop("'derivs' has to be a positive integer.")

    ## check and reformat 'degree'
    if ((degree <- as.integer(degree)) < 0)
        stop("'degree' must be a nonnegative integer.")

    ## original degree of freedom from definition
    df0 <- degree + length(knots) + 1L
    df <- df0 - as.integer(! intercept)

    ## take care of possible NA's in `x`
    nax <- is.na(x)
    if (all(nax))
        stop("'x' cannot be all NA's!")
    nas <- any(nax)

    ## attribute knots for output
    knotsAttr <- if (is.null(knots)) numeric(0L) else knots

    ## for derivs > degree
    if (derivs > degree) {
        dMat <- matrix(0, nrow = length(x), ncol = df)
        if (nas)
            dMat[nax, ] <- NA
        tmp <- list(degree = degree,
                    knots = knotsAttr,
                    Boundary.knots = Boundary.knots,
                    intercept = intercept,
                    x = x, derivs = derivs)
        attributes(dMat) <- c(attributes(dMat), tmp)
        class(dMat) <- c("matrix", "bSpline2", "deriv")
        return(dMat)
    }

    ## take care of possible NA's in `x`
    nax <- is.na(x)
    xx <- x
    ## remove NA's in x
    if (nas <- any(nax))
        xx <- x[! nax]

    ## check Boundary.knots specified by users
    if (! missing(Boundary.knots))
        Boundary.knots <- sort(Boundary.knots[seq_len(2)])

    dMat <- bSpline(xx, knots = knots, degree = degree - derivs,
                    intercept = TRUE, Boundary.knots = Boundary.knots)

    ## derivative matrix
    for (iter in seq_len(derivs)) {
        ## define knot sequence according to the bases being differentiated
        ord <- degree - derivs + iter + 1L
        aKnots <- sort(c(rep(Boundary.knots, ord), knots))
        denom <- diff(aKnots, lag = ord - 1L)
        facVec <- ifelse(abs(denom) < .Machine$double.eps,
                         0, (ord - 1L) / denom)
        dMat0 <- cbind(0, dMat, 0)
        dMat <- sapply(seq_len(df0 - derivs + iter), function(a)
        {
            idx <- a : (a + 1L)
            tmpMat <- dMat0[, idx, drop = FALSE]
            facVec[idx[1L]] * tmpMat[, 1L, drop = FALSE] -
                facVec[idx[2L]] * tmpMat[, 2L, drop = FALSE]
        })
    }

    ## take care of intercept
    if (! intercept)
        dMat <- dMat[, - 1L, drop = FALSE]

    ## keep NA's as is
    if (nas) {
        nmat <- matrix(NA, length(nax), ncol(dMat))
        nmat[! nax, ] <- dMat
        dMat <- nmat
    }

    ## add colnames for consistency with returns from splines::bs
    colnames(dMat) <- as.character(seq_len(df))

    ## on attributes
    tmp <- list(degree = degree,
                knots = knotsAttr,
                Boundary.knots = Boundary.knots,
                intercept = intercept,
                x = x, derivs = derivs)
    attributes(dMat) <- c(attributes(dMat), tmp)
    class(dMat) <- c("matrix", "dbs")
    ## return
    dMat
}
