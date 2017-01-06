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
##' \code{deriv} method that obtains derivative of given order of spline
##' bases. The function is designed for objects generated from this package.
##' The function extracts necessary information about the input spline basis
##' matrix from its attributes. So the function will not work if some
##' attribute is not availbale. For internal knots, the derivative is defined
##' to be the left derivative.
##'
##' @name deriv
##'
##' @param expr \code{spline2} objects generated from this package.
##' @param derivs A positive integer specifying the order of
##' derivatives. By default, it is \code{1L} for the first derivative.
##' @param ... Other arguments for further usage.
##'
##' @return A matrix of dimension \code{length(x)} by
##' \code{df = degree + length(knots)} (plus one if intercept is included).
##' Attributes that correspond to the arguments specified are returned
##' for usage for other function in this package.
##' @references
##' De Boor, Carl. (1978). \emph{A practical guide to splines}.
##' Vol. 27. New York: Springer-Verlag.
##' @examples
##' library(splines2)
##' x <- c(seq(0, 1, 0.1), NA)
##' knots <- c(0.3, 0.5, 0.6)
##' ## integal of B-splines
##' ibsMat <- ibs(x, knots = knots, degree = 3L, intercept = TRUE)
##' ## the corresponding B-splines integrated
##' bsMat <- bSpline(x, knots = knots, degree = 3L, intercept = TRUE)
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
deriv.bSpline2 <- function(expr, derivs = 1L, ...)
{
    ## extract splines information from attributes
    x <- attr(expr, "x")
    degree <- attr(expr, "degree")
    knots <- attr(expr, "knots")
    intercept <- attr(expr, "intercept")
    Boundary.knots <- attr(expr, "Boundary.knots")

    ## call function dbs
    dMat <- dbs(x = x, derivs = derivs, knots = knots, degree = degree,
                intercept = intercept, Boundary.knots = Boundary.knots, ...)

    ## prepare for output
    class(dMat) <- c("matrix", "bSpline2", "deriv")
    dMat
}


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
        class(out) <- c("matrix", "ibs", "deriv")
        return(out)
    }

    ## for higher order of derivative
    x <- attr(expr, "x")
    degree <- attr(expr, "degree")
    knots <- attr(expr, "knots")
    intercept <- attr(expr, "intercept")
    Boundary.knots <- attr(expr, "Boundary.knots")

    ## call function dbs
    dMat <- dbs(x = x, derivs = derivs - 1L, knots = knots, degree = degree,
                intercept = intercept, Boundary.knots = Boundary.knots, ...)

    ## prepare for output
    attr(dMat, "derivs") <- derivs
    class(dMat) <- c("matrix", "ibs", "deriv")
    dMat
}


##' @rdname deriv
##' @export
deriv.mSpline <- function(expr, derivs = 1L, ...)
{
    ## extract splines information from attributes
    x <- attr(expr, "x")
    degree <- attr(expr, "degree")
    knots <- attr(expr, "knots")
    intercept <- attr(expr, "intercept")
    Boundary.knots <- attr(expr, "Boundary.knots")

    ## call function dbs
    dMat <- mSpline(x = x, knots = knots, degree = degree,
                    intercept = intercept, Boundary.knots = Boundary.knots,
                    derivs = derivs, ...)

    ## prepare for output
    class(dMat) <- c("matrix", "mSpline", "deriv")
    dMat
}



##' @rdname deriv
##' @export
deriv.iSpline <- function(expr, derivs = 1L, ...)
{
    ## quick check on derivs
    derivs <- as.integer(derivs)
    if (derivs < 1L)
        stop("'derivs' has to be a positive integer.")

    ## extract existing result from attributes for the first derivative
    dMat <- if (derivs == 1L)
                attr(expr, "msMat")
            else
                ## for derivative of higher order
                deriv.mSpline(expr = expr, derivs = derivs - 1L, ...)

    ## prepare for output
    attr(dMat, "derivs") <- derivs
    class(dMat) <- c("matrix", "iSpline", "deriv")
    dMat
}


##' @rdname deriv
##' @export
deriv.cSpline <- function(expr, derivs = 1L, ...)
{
    derivs <- as.integer(derivs)
    if (derivs < 1L)
        stop("'derivs' has to be a positive integer.")

    dMat <- if (derivs == 1L) {
                attr(expr, "isMat")
            } else if (derivs == 2L) {
                attr(expr, "msMat")
            } else
                deriv.mSpline(expr = expr, derivs = derivs - 2L, ...)

    attr(dMat, "derivs") <- derivs
    class(dMat) <- c("matrix", "cSpline", "deriv")
    dMat
}
