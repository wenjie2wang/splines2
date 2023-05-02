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

##' M-Spline Basis for Polynomial Splines
##'
##' Generates the basis matrix of regular M-spline, periodic M-spline, and the
##' corresponding integrals and derivatives.
##'
##' This function contains an implementation of the closed-form M-spline basis
##' based on the recursion formula given by Ramsay (1988) or periodic M-spline
##' basis following the procedure producing periodic B-splines given in Piegl
##' and Tiller (1997).  For monotone regression, one can use I-splines (see
##' \code{\link{iSpline}}) instead of M-splines.  For shape-restricted
##' regression, see Wang and Yan (2021) for examples.
##'
##' The function \code{msp()} is an alias of to encourage the use in a model
##' formula.
##'
##' @inheritParams bSpline
##' @inherit bSpline return
##'
##' @references
##' Ramsay, J. O. (1988). Monotone regression splines in action.
##' \emph{Statistical science}, 3(4), 425--441.
##'
##' Piegl, L., & Tiller, W. (1997). \emph{The NURBS book}. Springer Science &
##' Business Media.
##'
##' Wang, W., & Yan, J. (2021). \emph{Shape-restricted regression splines with R
##' package splines2}. Journal of Data Science, 19(3),498--517.
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
    out <- rcpp_mSpline(
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
    class(out) <- c("MSpline", "splines2", "matrix")
    out
}

##' @rdname mSpline
##' @export
msp <- mSpline
