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
##' This function produces the derivative of given order of B-splines.  It is an
##' implementation of the close form derivative of B-spline basis based on
##' recursion relation.  At knots, the derivative is defined to be the right
##' derivative.
##'
##' The function is similar with \code{\link[splines]{splineDesign}}. However,
##' it provides a more user-friendly interface, a more considerate \code{NA}'s
##' handling.  Internally, it calls \code{\link{bSpline}} and generates a basis
##' matrix for representing the family of piecewise polynomials and their
##' corresponding derivative with the specified interior knots and degree,
##' evaluated at the values of \code{x}. The function \code{splineDesign} in
##' \code{splines} package can also be used to calculate derivative of
##' B-splines.
##'
##' @usage
##' dbs(x, derivs = 1L, df = NULL, knots = NULL, degree = 3L,
##'     intercept = FALSE, Boundary.knots = range(x, na.rm = TRUE), ...)
##'
##' @param x The predictor variable.  Missing values are allowed and will be
##'     kept and returned as they were.
##' @param derivs A positive integer specifying the order of derivative.  By
##'     default, it is \code{1L} for the first derivative.
##' @param df Degrees of freedom of the B-spline basis to be differentiated.
##'     One can specify \code{df} rather than \code{knots}, then the function
##'     chooses "df - degree" (minus one if there is an intercept) knots at
##'     suitable quantiles of \code{x} (which will ignore missing values).  The
##'     default, \code{NULL}, corresponds to no inner knots, i.e.,
##'     "degree - intercept".
##' @param knots The internal breakpoints that define the B-spline basis to be
##'     differentiated.  The default is \code{NULL}, which results in a basis
##'     for ordinary polynomial regression.  Typical values are the mean or
##'     median for one knot, quantiles for more knots.  See also
##'     \code{Boundary.knots}.
##' @param degree Non-negative integer degree of the piecewise polynomial to be
##'     differentiated. The default value is 3 for the integral of cubic
##'     B-splines.
##' @param intercept If \code{TRUE}, an intercept is included in the basis;
##'     Default is \code{FALSE}.
##' @param Boundary.knots Boundary points at which to anchor the B-spline basis
##'     to be differentiated. By default, they are the range of the
##'     non-\code{NA} data.  If both \code{knots} and \code{Boundary.knots} are
##'     supplied, the basis parameters do not depend on \code{x}.
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
##' x <- seq.int(0, 1, 0.01)
##' knots <- c(0.2, 0.4, 0.7)
##' ## the second derivative of cubic B-splines with three internal knots
##' dMat <- dbs(x, derivs = 2L, knots = knots, intercept = TRUE)
##'
##' ## compare with the results from splineDesign
##' ord <- attr(dMat, "degree") + 1L
##' bKnots <- attr(dMat, "Boundary.knots")
##' aKnots <- c(rep(bKnots[1L], ord), knots, rep(bKnots[2L], ord))
##' res <- splines::splineDesign(aKnots, x = x, derivs = 2L)
##' stopifnot(all.equal(res, dMat, check.attributes = FALSE))
##' @seealso
##' \code{\link{predict.dbs}} for evaluation at given (new) values;
##' \code{\link{deriv.dbs}} for derivative method;
##' \code{\link{bSpline}} for B-splines;
##' \code{\link{ibs}} for integral of B-splines.
##' @importFrom stats quantile
##' @export
dbs <- function(x, derivs = 1L, df = NULL, knots = NULL, degree = 3L,
                intercept = FALSE, Boundary.knots = range(x, na.rm = TRUE), ...)
{
    ## check order of derivative
    derivs <- as.integer(derivs)
    if (derivs < 1L)
        stop("'derivs' has to be a positive integer.")

    ## check and reformat 'degree'
    if ((degree <- as.integer(degree)) < 0)
        stop("'degree' must be a nonnegative integer.")

    ## sort and remove possible NA's in internal knots if exist
    if (length(knots)) {
        knots <- sort(knots)
        tmp <- is.na(knots)
        if (any(tmp)) {
            omit <- seq_along(knots)[tmp]
            knots <- knots[- omit]
        }
    }

    ## take care of possible NA's in `x`
    nax <- is.na(x)
    if (all(nax))
        stop("'x' cannot be all NA's!")
    nas <- any(nax)
    ## remove NA's
    xx <- if (nas) x[! nax] else x

    ## check Boundary.knots specified by users
    outside <- rep(FALSE, length(xx))
    if (! missing(Boundary.knots)) {
        Boundary.knots <- sort(Boundary.knots[seq_len(2)])
        outside <- (xx < Boundary.knots[1L]) | (xx > Boundary.knots[2L])
    }

    ## determine knots from df if missing
    inter <- as.integer(intercept)
    if (! is.null(df)) {
        df0 <- length(knots) + degree + inter
        if (tmp <- (df < df0))
            warning(gettextf("'df' was too small; have used %d",
                             df0), domain = NA)

        df <- ifelse(tmp, df0, df)
        nKnots <- df - degree - inter
        if (is.null(knots) && nKnots > 0) {
            quans <- seq.int(from = 0, to = 1,
                             length.out = nKnots + 2L)[- c(1L, nKnots + 2L)]
            knots <- stats::quantile(xx[! outside], quans)
        }
    }
    ## update degree of freedom from inputs
    df0 <- length(knots) + degree + 1L
    df <- df0 - 1L + inter

    ## attribute knots for output
    knotsAttr <- if (is.null(knots)) numeric(0L) else knots

    ## for derivs > degree
    if (derivs > degree) {
        ## df == 0, i.e., no basis returned
        if (! df)
            warning("Degree of freedom is zero.")
        dMat <- matrix(0, nrow = length(x), ncol = df)
        if (nas)
            dMat[nax, ] <- NA
        tmp <- list(degree = degree,
                    knots = knotsAttr,
                    Boundary.knots = Boundary.knots,
                    intercept = intercept,
                    x = x, derivs = derivs)
        attributes(dMat) <- c(attributes(dMat), tmp)
        class(dMat) <- c("matrix", "dbs")
        return(dMat)
    }

    ## B-spline bases
    dMat <- bSpline(xx, knots = knots, degree = degree - derivs,
                    intercept = TRUE, Boundary.knots = Boundary.knots, ...)

    ## derivative matrix
    for (iter in seq_len(derivs)) {
        ## define knot sequence according to the bases being differentiated
        ord <- degree - derivs + iter + 1L
        aKnots <- sort(c(rep(Boundary.knots, ord), knots))
        denom <- diff(aKnots, lag = ord - 1L)
        facVec <- ifelse(denom > 0, (ord - 1L) / denom, 0)
        dMat0 <- cbind(0, dMat, 0)
        dMat <- sapply(seq_len(df0 - derivs + iter), function(a)
        {
            idx <- a : (a + 1L)
            tmpMat <- dMat0[, idx, drop = FALSE]
            facVec[idx[1L]] * tmpMat[, 1L, drop = FALSE] -
                facVec[idx[2L]] * tmpMat[, 2L, drop = FALSE]
        })
    }

    ## recover dimension after sapply
    if (! is.matrix(dMat))
        dMat <- matrix(dMat, nrow = 1L)

    ## take care of intercept
    if (! intercept)
        dMat <- dMat[, - 1L, drop = FALSE]

    ## keep NA's as is
    if (nas) {
        nmat <- matrix(NA, length(nax), df)
        nmat[! nax, ] <- dMat
        dMat <- nmat
    }

    ## add dimnames for consistency with returns from splines::bs
    row.names(dMat) <- names(x)
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
