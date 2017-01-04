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


##' Derivative of Splines
##'
##' \code{deriv} method that obtains derivative of given number of spline
##' bases.  The function is designed for \code{splines2} objects generated
##' from this package.
##'
##' @name deriv
##'
##' @param expr \code{spline2} objects generated from this package.
##' @param derivs A positive integer specifying the number of
##' derivatives. By default, the
##' @param ... Other arguments for further usage.
##'
##' @return A matrix of dimension \code{length(x)} by
##' \code{df = degree + length(knots)} (plus one if intercept is included).
##' Attributes that correspond to the arguments specified are returned
##' for usage for other function in this package.
##' @examples
##' library(splines2)
##' x <- c(seq(0, 1, 0.1), NA)
##' knots <- c(0.3, 0.5, 0.6)               # set up knots
##' ## integal of B-splines
##' ibsMat <- ibs(x, knots = knots, degree = 3L, intercept = TRUE)
##' ## the corresponding B-splines integrated
##' bsMat <- ibs(x, knots = knots, degree = 3L, intercept = TRUE)
##' ## the first derivative should equal bsMat
##' d1Mat <- deriv(ibsMat)
##' all.equal(bsMat, d1Mat, check.attributes = FALSE)
##' ## the second derivative shoud equal the first derivative of bsMat
##' dbsMat <- deriv(bsMat)
##' d2Mat <- deriv(ibsMat, derivs = 2L)
##' all.equal(dbsMat, d2Mat, check.attributes = FALSE)
##' @seealso
##' \code{\link{bSpline}} for B-splines;
##' \code{\link{ibs}} for integral of B-spline basis;
##' \code{\link{mSpline}} for M-spline basis;
##' \code{\link{iSpline}} for I-spline basis;
##' \code{\link{cSpline}} for C-spline basis.
##' @importFrom stats deriv
NULL


##' @rdname deriv
##' @export
deriv.ibs <- function(expr, derivs = 1L, ...)
{
    ## quick check on derivs
    derivs <- as.integer(derivs)
    if (derivs < 1L)
        stop("'derivs' has to be a positive integer.")

    ## if first derivative, take result from existing attribute
    if (derivs == 1L) {
        out <- attr(expr, "bsMat")
        attr(out, "derivs") <- derivs
        class(out) <- c("deriv", "basis", "matrix")
        return(out)
    }

    ## otherwise
    newx <- attr(expr, "x")
    bsMat <- predict.bSpline2(expr, newx = newx)
    deriv.bSpline2(bsMat, derivs = derivs - 1L, ...)
}


##' @rdname deriv
##' @export
deriv.iSpline <- function(expr, derivs = 1L, ...)
{
    derivs <- as.integer(derivs)
    if (derivs < 1L)
        stop("'derivs' has to be a positive integer.")

    ## take care of possible NA's in `x`
    x <- attr(expr, "x")
    nax <- is.na(x)
    nas <- any(nax)

    degree <- attr(expr, "degree")
    if (derivs > degree + 1L) {
        expr[, ] <- 0
        if (nas)
            expr[nax, ] <- NA
        attr(expr, "derivs") <- derivs
        class(expr) <- c("deriv", "basis", "matrix")
        return(expr)
    }

    if (derivs == 1L)
        return(attr(expr, "msMat"))

    stop("under development")
}


##' @rdname deriv
##' @export
deriv.cSpline <- function(expr, derivs = 1L, ...)
{
    derivs <- as.integer(derivs)
    if (derivs < 1L)
        stop("'derivs' has to be a positive integer.")

    ## take care of possible NA's in `x`
    x <- attr(expr, "x")
    nax <- is.na(x)
    nas <- any(nax)

    degree <- attr(expr, "degree")
    if (derivs > degree + 1L) {
        expr[, ] <- 0
        if (nas)
            expr[nax, ] <- NA
        attr(expr, "derivs") <- derivs
        class(expr) <- c("deriv", "basis", "matrix")
        return(expr)
    }

    if (derivs == 1L)
        return(attr(expr, "isMat"))
    if (derivs == 2L)
        return(attr(expr, "msMat"))

    stop("under development")
}


##' @rdname deriv
##' @export
deriv.bSpline2 <- function(expr, derivs = 1L, ...)
{
    derivs <- as.integer(derivs)
    if (derivs < 1L)
        stop("'derivs' has to be a positive integer.")

    ## extract splines information from attributes
    degree <- attr(expr, "degree")
    knots <- attr(expr, "knots")
    intercept <- attr(expr, "intercept")
    df <- degree + length(knots) + as.integer(intercept)
    x <- attr(expr, "x")

    ## take care of possible NA's in `x`
    xx <- x
    nax <- is.na(x)
    ## remove NA's in x
    if (nas <- any(nax))
        xx <- x[! nax]

    ## for derivs > degree
    if (derivs > degree) {
        expr[, ] <- 0
        if (nas)
            expr[nax, ] <- NA
        attr(expr, "derivs") <- derivs
        class(expr) <- c("deriv", "basis", "matrix")
        return(expr)
    }

    ## remove possible NA's
    dMat <- expr
    attr(dMat, "x") <- xx

    ## call internal function derivBs
    dMat <- derivBs(dMat, derivs = derivs, ...)

    ## keep NA's as is
    if (nas) {
        nmat <- matrix(NA, length(nax), ncol(dMat))
        nmat[! nax, ] <- dMat
        dMat <- nmat
    }

    ## add colnames for consistency with returns from splines::bs
    colnames(dMat) <- as.character(seq_len(df))

    ## on attributes
    attributes(dMat) <- c(attributes(expr), list(derivs = derivs))
    class(dMat) <- c("deriv", "basis", "matrix")

    ## return
    dMat
}


### internal function ==========================================================
## function computing the derivative of B-spline bases
derivBs <- function(bsMat, derivs = 1L, ...)
{
    ## function for internally usage only
    ## checking procedure is omitted for performance

    ## extract splines information from attributes
    degree <- attr(bsMat, "degree")
    knots <- attr(bsMat, "knots")
    intercept <- attr(bsMat, "intercept")
    bKnots <- attr(bsMat, "Boundary.knots")

    ## original degree of freedom from definition
    df0 <- degree + length(knots) + 1L

    ## linear bases
    x <- attr(bsMat, "x")
    dMat <- bSpline(x, knots = knots, degree = degree - derivs,
                    intercept = TRUE, Boundary.knots = bKnots)

    ## derivative matrix
    for (iter in seq_len(derivs)) {
        ## define knot sequence according to the bases being differentiated
        ord <- degree - derivs + iter + 1L
        aKnots <- sort(c(rep(bKnots, ord), knots))
        denom1 <- diff(aKnots, lag = ord - 1L)
        facVec <- ifelse(denom1 > 0, (ord - 1L) / denom1, 0)
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

    ## return
    dMat
}
