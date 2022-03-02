##
## R package splines2 by Wenjie Wang and Jun Yan
## Copyright (C) 2016-2022
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

##' M-Spline Basis for Polynomial Splines
##'
##' Generates the basis matrix of regular M-spline, periodic M-spline, and the
##' corresponding integrals and derivatives.
##'
##' This function contains an implementation of the closed-form M-spline basis
##' based on the recursion formula given by Ramsay (1988) or periodic M-spline
##' basis following the procedure producing periodic B-splines given in Piegl
##' and Tiller (1997).  For monotone regression, one can use I-splines (see
##' \code{\link{iSpline}}) instead of M-splines.
##'
##' @inheritParams bSpline
##'
##' @param df Degree of freedom that equals to the column number of the returned
##'     matrix.  One can specify \code{df} rather than \code{knots}.  For
##'     M-splines, the function chooses \code{df - degree -
##'     as.integer(intercept)} internal knots at suitable quantiles of \code{x}
##'     ignoring missing values and those \code{x} outside of the boundary.  For
##'     periodic M-spline (\code{periodic = TRUE}), \code{df -
##'     as.integer(intercept)} internal knots will be chosen at suitable
##'     quantiles of \code{x} relative to the beginning of the cyclic intervals
##'     they belong to (see Examples) and the number of internal knots must be
##'     greater or equal to the specified \code{degree - 1}.  If internal knots
##'     are specified via \code{knots}, the specified \code{df} will be ignored.
##' @param knots The internal breakpoints that define the splines.  The default
##'     is \code{NULL}, which results in a basis for ordinary polynomial
##'     regression.  Typical values are the mean or median for one knot,
##'     quantiles for more knots.  For periodic splines (\code{periodic =
##'     TRUE}), the number of knots must be greater or equal to the specified
##'     \code{degree - 1}.
##' @param Boundary.knots Boundary points at which to anchor the splines.  By
##'     default, they are the range of \code{x} excluding \code{NA}.  If both
##'     \code{knots} and \code{Boundary.knots} are supplied, the basis
##'     parameters do not depend on \code{x}. Data can extend beyond
##'     \code{Boundary.knots}.  For periodic splines (\code{periodic = TRUE}),
##'     the specified boundary knots define the cyclic interval.
##' @param periodic A logical value.  If \code{TRUE}, the periodic splines will
##'     be returned instead of regular M-splines.  The default value is
##'     \code{FALSE}.
##' @param derivs A nonnegative integer specifying the order of derivatives of
##'     M-splines. The default value is \code{0L} for M-spline basis functions.
##' @param integral A logical value.  If \code{TRUE}, the corresponding
##'     integrals of spline basis functions will be returned.  The default value
##'     is \code{FALSE}.  For periodic splines, the integral of each basis is
##'     integrated from the left boundary knot.
##'
##' @return A numeric matrix of \code{length(x)} rows and \code{df} columns if
##'     \code{df} is specified.  If \code{knots} are specified instead, the
##'     output matrix will consist of \code{length(knots) + degree +
##'     as.integer(intercept)} columns if \code{periodic = FALSE}, or
##'     \code{length(knots) + as.integer(intercept)} columns if \code{periodic =
##'     TRUE}.  Attributes that correspond to the arguments specified are
##'     returned for usage of other functions in this package.
##'
##' @references
##' Ramsay, J. O. (1988). Monotone regression splines in action.
##' \emph{Statistical science}, 3(4), 425--441.
##'
##' Piegl, L., & Tiller, W. (1997). \emph{The NURBS book}. Springer Science \&
##' Business Media.
##'
##' @example inst/examples/ex-mSpline.R
##'
##' @seealso
##' \code{\link{bSpline}} for B-splines;
##' \code{\link{iSpline}} for I-splines;
##' \code{\link{cSpline}} for C-splines.
##'
##' @export
mSpline <- function(x, df = NULL, knots = NULL, degree = 3L,
                    intercept = FALSE, Boundary.knots = NULL,
                    periodic = FALSE, derivs = 0L, integral = FALSE,
                    ...)
{
    ## check inputs
    if ((derivs <- as.integer(derivs)) < 0) {
        stop("The 'derivs' must be a nonnegative integer.")
    }
    if ((degree <- as.integer(degree)) < 0)
        stop("The 'degree' must be a nonnegative integer.")
    if (is.null(df)) {
        df <- 0L
    } else {
        df <- as.integer(df)
        if (df < 0) {
            stop("The 'df' must be a nonnegative integer.")
        } else if (periodic && df < degree) {
            stop("The 'df' must be >= 'degree' for periodic spline basis.")
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
    out <- if (periodic) {
               rcpp_periodic_mSpline(
                   x = xx,
                   df = df,
                   degree = degree,
                   internal_knots = knots,
                   boundary_knots = Boundary.knots,
                   derivs = derivs,
                   integral = integral,
                   complete_basis = intercept
               )
           } else {
               rcpp_mSpline(
                   x = xx,
                   df = df,
                   degree = degree,
                   internal_knots = knots,
                   boundary_knots = Boundary.knots,
                   derivs = derivs,
                   integral = integral,
                   complete_basis = intercept
               )
           }
    ## throw warning if any x is outside of the boundary
    b_knots <- attr(out, "Boundary.knots")
    if (! periodic && any((xx < b_knots[1L]) | (xx > b_knots[2L]))) {
        warning(wrapMessages(
            "Some 'x' values beyond boundary knots",
            "may cause ill-conditioned basis functions."
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
    class(out) <- c("matrix", "mSpline")
    out
}
