##
## R package splines2 by Wenjie Wang and Jun Yan
## Copyright (C) 2016-2023
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

##' B-Spline Basis for Polynomial Splines
##'
##' Generates the B-spline basis matrix representing the family of piecewise
##' polynomials with the specified interior knots, degree, and boundary knots,
##' evaluated at the values of \code{x}.
##'
##' This function extends the \code{bs()} function in the \code{splines} package
##' for B-spline basis by allowing piecewise constant (left-closed and
##' right-open except on the right boundary) spline basis of degree zero.
##'
##' @param x The predictor variable.  Missing values are allowed and will be
##'     returned as they are.
##' @param df Degree of freedom that equals to the column number of the returned
##'     matrix.  One can specify \code{df} rather than \code{knots}, then the
##'     function chooses \code{df - degree - as.integer(intercept)} internal
##'     knots at suitable quantiles of \code{x} ignoring missing values and
##'     those \code{x} outside of the boundary.  If internal knots are specified
##'     via \code{knots}, the specified \code{df} will be ignored.
##' @param knots The internal breakpoints that define the splines.  The default
##'     is \code{NULL}, which results in a basis for ordinary polynomial
##'     regression.  Typical values are the mean or median for one knot,
##'     quantiles for more knots.
##' @param degree A nonnegative integer specifying the degree of the piecewise
##'     polynomial. The default value is \code{3} for cubic splines. Zero degree
##'     is allowed for piecewise constant basis functions.
##' @param intercept If \code{TRUE}, the complete basis matrix will be returned.
##'     Otherwise, the first basis will be excluded from the output.
##' @param Boundary.knots Boundary points at which to anchor the splines.  By
##'     default, they are the range of \code{x} excluding \code{NA}.  If both
##'     \code{knots} and \code{Boundary.knots} are supplied, the basis
##'     parameters do not depend on \code{x}. Data can extend beyond
##'     \code{Boundary.knots}.
##' @param derivs A nonnegative integer specifying the order of derivatives of
##'     B-splines. The default value is \code{0L} for B-spline basis functions.
##' @param integral A logical value.  If \code{TRUE}, the corresponding
##'     integrals of spline basis functions will be returned.  The default value
##'     is \code{FALSE}.
##' @param ... Optional arguments that are not used.
##'
##' @return A numeric matrix of \code{length(x)} rows and \code{df} columns if
##'     \code{df} is specified or \code{length(knots) + degree +
##'     as.integer(intercept)} columns if \code{knots} are specified instead.
##'     Attributes that correspond to the arguments specified are returned
##'     mainly for other functions in this package.
##'
##' @references
##' De Boor, Carl. (1978). \emph{A practical guide to splines}.
##' Vol. 27. New York: Springer-Verlag.
##'
##' @example inst/examples/ex-bSpline.R
##'
##' @seealso
##' \code{\link{dbs}} for derivatives of B-splines;
##' \code{\link{ibs}} for integrals of B-splines;
##'
##' @export
bSpline <- function(x, df = NULL, knots = NULL, degree = 3L,
                    intercept = FALSE, Boundary.knots = NULL,
                    derivs = 0L, integral = FALSE, ...)
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
    out <- rcpp_bSpline(
        x = xx,
        df = df,
        degree = degree,
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
            "may cause ill-conditioned basis functions."
        ))
    }
    ## keep NA's as is
    if (nas) {
        nmat <- matrix(NA_real_, length(nax), ncol(out))
        nmat[! nax, ] <- out
        saved_attr <- attributes(out)
        saved_attr$dim[1] <- length(nax)
        out <- nmat
        attributes(out) <- saved_attr
        attr(out, "x") <- x
    }
    ## add dimnames for consistency with bs returns
    name_x <- names(x)
    if (! is.null(name_x)) {
        row.names(out) <- name_x
    }
    ## add class
    class(out) <- c("matrix", "bSpline2", "splines2")
    ## return
    out
}
