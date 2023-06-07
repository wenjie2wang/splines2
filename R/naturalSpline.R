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

##' Natural Cubic Spline Basis for Polynomial Splines
##'
##' Functions \code{naturalSpline()} and \code{nsk()} generate the natural cubic
##' spline basis functions, the corresponding derivatives or integrals (from the
##' left boundary knot).  Both of them are different from \code{splines::ns()}.
##' However, for a given model fitting procedure, using different variants of
##' spline basis functions should result in identical prediction values.  The
##' coefficient estimates of the spline basis functions returned by \code{nsk()}
##' are more interpretable compared to \code{naturalSpline()} or
##' \code{splines::ns()} .
##'
##' The constructed spline basis functions from \code{naturalSpline()} are
##' nonnegative within boundary with the second derivatives being zeros at
##' boundary knots.  The implementation utilizes the close-form null space that
##' can be derived from the recursive formula for the second derivatives of
##' B-splines.  The function \code{nsp()} is an alias of \code{naturalSpline()}
##' to encourage the use in a model formula.
##'
##' The function \code{nsk()} produces another variant of natural cubic spline
##' matrix such that only one of the basis functions is nonzero and takes a
##' value of one at every boundary and internal knot.  As a result, the
##' coefficients of the resulting fit are the values of the spline function at
##' the knots, which makes it easy to interpret the coefficient estimates.  In
##' other words, the coefficients of a linear model will be the heights of the
##' function at the knots if \code{intercept = TRUE}.  If \code{intercept =
##' FALSE}, the coefficients will be the change in function value between each
##' knot.  This implementation closely follows the function \code{nsk()} of the
##' \pkg{survival} package (version 3.2-8).  The idea corresponds directly to
##' the physical implementation of a spline by passing a flexible strip of wood
##' or metal through a set of fixed points, which is a traditional way to create
##' smooth shapes for things like a ship hull.
##'
##' The returned basis matrix can be obtained by transforming the corresponding
##' B-spline basis matrix with the matrix \code{H} provided in the attribute of
##' the returned object.  Each basis is assumed to follow a linear trend for
##' \code{x} outside of boundary.  A similar implementation is provided by
##' \code{splines::ns}, which uses QR decomposition to find the null space of
##' the second derivatives of B-spline basis at boundary knots.  See
##' Supplementray Materials of Wang and Yan (2021) for a more detailed
##' introduction.
##'
##' @name naturalSpline
##'
##' @inheritParams bSpline
##'
##' @param df Degree of freedom that equals to the column number of returned
##'     matrix.  One can specify \code{df} rather than \code{knots}, then the
##'     function chooses \code{df - 1 - as.integer(intercept)} internal knots at
##'     suitable quantiles of \code{x} ignoring missing values and those
##'     \code{x} outside of the boundary.  Thus, \code{df} must be greater than
##'     or equal to \code{2}.  If internal knots are specified via \code{knots},
##'     the specified \code{df} will be ignored.
##' @param derivs A nonnegative integer specifying the order of derivatives of
##'     natural splines. The default value is \code{0L} for the spline basis
##'     functions.
##' @param integral A logical value.  The default value is \code{FALSE}.  If
##'     \code{TRUE}, this function will return the integrated natural splines
##'     from the left boundary knot.
##' @param trim The fraction (0 to 0.5) of observations to be trimmed from each
##'     end of \code{x} before placing the default internal and boundary knots.
##'     This argument will be ignored if \code{Boundary.knots} is specified.
##'     The default value is \code{0} for backward compatability, which sets the
##'     boudary knots as the range of |code{x}.  If a positive fraction is
##'     specified, the default boundary knots will be equivalent to
##'     \code{quantile(x, probs = c(trim, 1 - trim), na.rm = TRUE)}, which can
##'     be a more sensible choice in practice due to the existence of outliers.
##'     The default internal knots are placed within the boundary afterwards.
##'
##' @return A numeric matrix of \code{length(x)} rows and \code{df}
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
NULL

## engine function
.engine_nsp <- function(x, df = NULL, knots = NULL,
                        intercept = FALSE, Boundary.knots = NULL,
                        trim = 0, derivs = 0L, integral = FALSE, ...,
                        .FUN = c("nsp", "nsk"))
{
    ## check inputs
    if ((derivs <- as.integer(derivs)) < 0) {
        stop("The 'derivs' must be a nonnegative integer.")
    }
    if (is.null(df)) {
        df <- 0L
    } else {
        df <- as.integer(df)
        if (df < 2) {
            stop("The 'df' must be >= 2.")
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
    .FUN <- match.arg(.FUN, choices = c("nsp", "nsk"))
    if (.FUN == "nsp") {
        out <- rcpp_naturalSpline(
            x = xx,
            df = df,
            internal_knots = knots,
            boundary_knots = Boundary.knots,
            trim = trim,
            complete_basis = intercept,
            derivs = derivs,
            integral = integral
        )
        nsp_class <- "NaturalSpline"
    } else {
        out <- rcpp_nsk(
            x = xx,
            df = df,
            internal_knots = knots,
            boundary_knots = Boundary.knots,
            trim = trim,
            complete_basis = intercept,
            derivs = derivs,
            integral = integral
        )
        nsp_class <- "NaturalSplineK"
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
    class(out) <- c(nsp_class, "splines2", "matrix")
    out
}

##' @rdname naturalSpline
##' @export
naturalSpline <- function(x, df = NULL, knots = NULL,
                          intercept = FALSE, Boundary.knots = NULL,
                          trim = 0, derivs = 0L, integral = FALSE, ...)
{
    .engine_nsp(
        x = x,
        df = df,
        knots = knots,
        intercept = intercept,
        Boundary.knots = Boundary.knots,
        trim = trim,
        derivs = derivs,
        integral = integral,
        .FUN = "nsp"
    )
}


##' @rdname naturalSpline
##' @export
nsp <- naturalSpline


##' @rdname naturalSpline
##' @export
nsk <- function(x, df = NULL, knots = NULL,
                intercept = FALSE, Boundary.knots = NULL,
                trim = 0, derivs = 0L, integral = FALSE, ...)
{
    .engine_nsp(
        x = x,
        df = df,
        knots = knots,
        intercept = intercept,
        Boundary.knots = Boundary.knots,
        trim = trim,
        derivs = derivs,
        integral = integral,
        .FUN = "nsk"
    )
}
