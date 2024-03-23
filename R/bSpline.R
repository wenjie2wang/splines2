##
## R package splines2 by Wenjie Wang and Jun Yan
## Copyright (C) 2016-2024
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
##' Generates the spline basis matrix for B-splines representing the family of
##' piecewise polynomials with the specified interior knots, degree, and
##' boundary knots, evaluated at the values of \code{x}.
##'
##' This function extends the \code{bs()} function in the \code{splines} package
##' for B-spline basis functions by allowing piecewise constant (left-closed and
##' right-open except on the right boundary) spline basis of degree zero.  In
##' addition, the function provides derivatives or integrals of the B-spline
##' basis functions when one specifies the arguments \code{derivs} or
##' \code{integral} appropriately.  The function constructs periodic B-splines
##' when \code{periodic} is \code{TRUE}.  All the implementations are based on
##' the closed-form recursion formula following De Boor (1978) and Wang and Yan
##' (2021).
##'
##' The functions \code{ibs()} and \code{dbs()} are provided for convenience.
##' The former provides the integrals of B-splines and is equivalent to
##' \code{bSpline()} with \code{integral = TRUE}.  The latter produces the
##' derivatives of given order of B-splines and is equivalent to
##' \code{bSpline()} with default \code{derivs = 1}.  The function \code{bsp()}
##' is an alias of to encourage the use in a model formula.
##'
##' @param x The predictor variable.  Missing values are allowed and will be
##'     returned as they are.
##' @param df Degree of freedom that equals to the column number of the returned
##'     matrix.  One can specify \code{df} rather than \code{knots}, then the
##'     function chooses \code{df - degree - as.integer(intercept)} internal
##'     knots at suitable quantiles of \code{x} ignoring missing values and
##'     those \code{x} outside of the boundary.  For periodic splines, \code{df
##'     - as.integer(intercept)} internal knots will be chosen at suitable
##'     quantiles of \code{x} relative to the beginning of the cyclic intervals
##'     they belong to (see Examples) and the number of internal knots must be
##'     greater or equal to the specified \code{degree - 1}.  If internal knots
##'     are specified via \code{knots}, the specified \code{df} will be ignored.
##' @param knots The internal breakpoints that define the splines.  The default
##'     is \code{NULL}, which results in a basis for ordinary polynomial
##'     regression.  Typical values are the mean or median for one knot,
##'     quantiles for more knots.  For periodic splines, the number of knots
##'     must be greater or equal to the specified \code{degree - 1}.
##'     Duplicated internal knots are not allowed.
##' @param degree A nonnegative integer specifying the degree of the piecewise
##'     polynomial. The default value is \code{3} for cubic splines. Zero degree
##'     is allowed for piecewise constant basis functions.
##' @param intercept If \code{TRUE}, the complete basis matrix will be returned.
##'     Otherwise, the first basis will be excluded from the output.
##' @param Boundary.knots Boundary points at which to anchor the splines.  By
##'     default, they are the range of \code{x} excluding \code{NA}.  If both
##'     \code{knots} and \code{Boundary.knots} are supplied, the basis
##'     parameters do not depend on \code{x}. Data can extend beyond
##'     \code{Boundary.knots}.  For periodic splines, the specified bounary
##'     knots define the cyclic interval.
##' @param periodic A logical value.  If \code{TRUE}, the periodic splines will
##'     be returned.  The default value is \code{FALSE}.
##' @param derivs A nonnegative integer specifying the order of derivatives of
##'     splines basis function.  The default value is \code{0}.
##' @param integral A logical value.  If \code{TRUE}, the corresponding
##'     integrals of spline basis functions will be returned.  The default value
##'     is \code{FALSE}. For periodic splines, the integral of each basis is
##'     integrated from the left boundary knot.
##' @param warn.outside A logical value indicating if a warning should be thrown
##'     out when any \code{x} is outside the boundary.  This option can also be
##'     set through \code{options("splines2.warn.outside")} after the package is
##'     loaded.
##' @param ... Optional arguments that are not used.
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
##' De Boor, Carl. (1978). \emph{A practical guide to splines}.
##' Vol. 27. New York: Springer-Verlag.
##'
##' Wang, W., & Yan, J. (2021). \emph{Shape-restricted regression splines with R
##' package splines2}. Journal of Data Science, 19(3),498--517.
##'
##' @example inst/examples/ex-bSpline.R
##'
##' @seealso
##' \code{\link{knots}} for extracting internal and boundary knots.
##'
##' @export
bSpline <- function(x, df = NULL, knots = NULL, degree = 3L,
                    intercept = FALSE, Boundary.knots = NULL,
                    periodic = FALSE, derivs = 0L, integral = FALSE,
                    warn.outside = getOption("splines2.warn.outside", TRUE),
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
        complete_basis = intercept,
        periodic = periodic,
        derivs = derivs,
        integral = integral
    )
    ## throw warning if any x is outside of the boundary
    b_knots <- attr(out, "Boundary.knots")
    if (warn.outside && ! periodic &&
        any((xx < b_knots[1L]) | (xx > b_knots[2L]))) {
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
    class(out) <- c("BSpline", "splines2", "matrix")
    ## return
    out
}

##' @rdname bSpline
##' @export
ibs <- function(x, df = NULL, knots = NULL, degree = 3,
                intercept = FALSE, Boundary.knots = NULL, ...)
{
    bSpline(x = x, df = df, knots = knots, degree = degree,
            intercept = intercept, Boundary.knots = Boundary.knots,
            integral = TRUE, ...)
}

##' @rdname bSpline
##' @export
dbs <- function(x, derivs = 1L, df = NULL, knots = NULL, degree = 3,
                intercept = FALSE, Boundary.knots = NULL, ...)
{
    bSpline(x = x, df = df, knots = knots, degree = degree,
            intercept = intercept, Boundary.knots = Boundary.knots,
            derivs = derivs, ...)
}

##' @rdname bSpline
##' @export
bsp <- bSpline
