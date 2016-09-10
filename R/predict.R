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


##' Evaluate a Spline Basis
##'
##' Evaluate a predefined spline basis at given values.
##'
##' These are methods for the generic function \code{predict} for objects
##' inheriting from class \code{\link{ibs}}, \code{\link{mSpline}}, or
##' \code{\link{iSpline}}. If \code{newx} is not given, the function returns
##' the input object.
##'
##' @name predict
##' @param object Objects of class \code{\link{ibs}}, \code{\link{mSpline}}, or
##' \code{\link{iSpline}} having attributes describing \code{knots},
##' \code{degree}, etc.
##' @param newx The \code{x} values at which evaluations are required.
##' @param ... Optional argument for future usage.
##'
##' @return An object just like the \code{object} input, except evaluated at
##' the new values of \code{x}.
##'
##' @examples
##' x <- seq(0, 1, 0.2)
##' knots <- c(0.3, 0.5, 0.6)
##' newX <- seq(0.1, 0.9, 0.2)
##'
##' ## for integral of B-spline
##' ibsMat <- ibs(x, knots = knots, degree = 2)
##' predict(ibsMat, newX)
##'
##' ## for M-spline
##' msMat <- mSpline(x, knots = knots, degree = 2)
##' predict(msMat, newX)
##'
##' ## for I-spline
##' imsMat <- iSpline(x, knots = knots, degree = 2)
##' predict(imsMat, newX)
##'
##' @seealso
##' \code{\link{ibs}} for integral of B-spline basis;
##' \code{\link{mSpline}} for M-spline basis;
##' \code{\link{iSpline}} for I-spline basis.
##'
##' @importFrom stats predict
NULL


##' @rdname predict
##' @export
predict.ibs <- function(object, newx, ...) {
    if (missing(newx)) return(object)
    a <- c(list(x = newx),
          attributes(object)[c("degree", "knots", "Boundary.knots",
                               "intercept")])
    do.call("ibs", a)
}


##' @rdname predict
##' @export
predict.mSpline <- function(object, newx, ...) {
    if (missing(newx)) return(object)
    a <- c(list(x = newx),
          attributes(object)[c("degree", "knots", "Boundary.knots",
                               "intercept")])
    do.call("mSpline", a)
}


##' @rdname predict
##' @export
predict.iSpline <- function(object, newx, ...) {
    if (missing(newx)) return(object)
    a <- c(list(x = newx),
          attributes(object)[c("degree", "knots", "Boundary.knots",
                               "intercept")])
    do.call("iSpline", a)
}


##' @rdname predict
##' @export
predict.cSpline <- function(object, newx, ...) {
    if (missing(newx)) return(object)
    a <- c(list(x = newx),
          attributes(object)[c("degree", "knots", "Boundary.knots",
                               "intercept")])
    do.call("cSpline", a)
}
