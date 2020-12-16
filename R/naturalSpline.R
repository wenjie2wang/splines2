##
## R package splines2 by Wenjie Wang and Jun Yan
## Copyright (C) 2016-2020
##
## This file is part of the R package splines2.
##
## The R package splines2 is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package splines2 is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##

##' Cubic Natural Spline Basis for Polynomial Splines
##'
##' Generates the nonnegative cubic natural spline basis matrix or the
##' corresponding derivatives of given order.
##'
##' It is an implementation of the natural spline basis based on B-spline basis.
##' The constructed spline bases are intended to be non-negative with second
##' derivatives being zeros at boundary knots.
##'
##' A similar implementation is provided by \code{splines::ns}, which uses QR
##' decomposition to find the null space of the second derivatives of B-spline
##' basis at boundary knots.  However, there is no guarantee that the resulting
##' bases are nonnegative over their support.
##'
##' @inheritParams bSpline
##'
##' @param derivs A non-negative integer specifying the order of derivatives of
##'     natural splines. The default value is \code{0L} for the spline bases.
##' @param integral A logical value.  If \code{TRUE}, the integrals of the
##'     natural spline bases will be returned.  The default value is
##'     \code{FALSE}.
##'
##' @return A numeric matrix with \code{length(x)} rows and \code{df}
##'     columns if \code{df} is specified or \code{length(knots) + 1 +
##'     as.integer(intercept)} columns if \code{knots} are specified instead.
##'     Attributes that correspond to the arguments specified are returned for
##'     usage of other functions in this package.
##'
##' @example inst/examples/ex-naturalSpline.R
##'
##' @seealso
##' \code{\link{bSpline}} for B-splines;
##' \code{\link{mSpline}} for M-splines;
##' \code{\link{iSpline}} for I-splines.
##'
##' @export
naturalSpline <- function(x, df = NULL, knots = NULL,
                          intercept = FALSE, Boundary.knots = NULL,
                          derivs = 0L, integral = FALSE,
                          ...)
{
    ## check inputs
    if ((derivs <- as.integer(derivs)) < 0) {
        stop("'derivs' must be a nonnegative integer.")
    }
    if (is.null(df)) {
        df <- 0L
    } else {
        df <- as.integer(df)
        if (df < 0) {
            stop("'df' must be a nonnegative integer.")
        }
    }
    knots <- null2num0(knots)
    Boundary.knots <- null2num0(Boundary.knots)
    ## take care of possible NA's in `x`
    nax <- is.na(x)
    if (all(nax)) {
        stop("The 'x' cannot be all NA's!")
    }
    ## remove NA's in x
    xx <- if (nas <- any(nax)) {
              x[! nax]
          } else {
              x
          }
    ## call the engine function
    out <- rcpp_naturalSpline(
        x = xx,
        df = df,
        internal_knots = knots,
        boundary_knots = Boundary.knots,
        derivs = derivs,
        integral = integral,
        complete_basis = intercept
    )
    ## throw warning if any x is outside of the boundary
    b_knots <- attr(out, "Boundary.knots")
    if (any((xx < b_knots[1L]) | (xx > b_knots[2L]))) {
        warning(wrapMessages(
            "Some 'x' values beyond boundary knots",
            "may cause ill-conditioned bases."
        ))
    }
    ## keep NA's as is
    if (nas) {
        nmat <- matrix(NA, length(nax), ncol(out))
        nmat[! nax, ] <- out
        saved_attr <- attributes(out)
        saved_attr$dim[1] <- length(nax)
        out <- nmat
        attributes(out) <- saved_attr
        attr(out, "x") <- x
    }
    ## add dimnames for consistency
    name_x <- names(x)
    if (! is.null(name_x)) {
        row.names(out) <- name_x
    }
    ## add class
    class(out) <- c("matrix", "naturalSpline")
    out
}