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


##' B-Spline Basis for Polynomial Splines
##'
##' This function generates the B-spline basis matrix for a polynomial spline.
##'
##' It is an augmented function of \code{\link[splines]{bs}} in package
##' \code{splines} for B-spline basis that allows piecewise constant spline
##' basis with zero degree. When the argument \code{degree} is greater than
##' zero, it internally calls \code{\link[splines]{bs}} and generates a basis
##' matrix for representing the family of piecewise polynomials with the
##' specified interior knots and degree, evaluated at the values of \code{x}.
##' The function has the same arguments with \code{\link[splines]{bs}} for ease
##' usage.
##'
##' @usage bSpline(x, df = NULL, knots = NULL, degree = 3, intercept = FALSE,
##'         Boundary.knots = range(x), ...)
##' @param x The predictor variable.  Missing values are allowed and will be
##' returned as they were.
##' @param df Degrees of freedom.  One can specify \code{df} rather than
##' \code{knots}, then the function chooses "df - degree"
##' (minus one if there is an intercept) knots at suitable quantiles of \code{x}
##' (which will ignore missing values).  The default, \code{NULL}, corresponds
##' to no inner knots, i.e., "degree - intercept".
##' @param knots The internal breakpoints that define the spline.  The
##' default is \code{NULL}, which results in a basis for ordinary
##' polynomial regression.  Typical values are the mean or median
##' for one knot, quantiles for more knots.  See also
##' \code{Boundary.knots}.
##' @param degree Non-negative integer degree of the piecewise polynomial. The
##' default value is 3 for cubic splines. Zero degree is allowed for this
##' function, which is the only difference compared with
##' \code{\link[splines]{bs}} in package \code{splines}.
##' @param intercept If \code{TRUE}, an intercept is included in the basis;
##' Default is \code{FALSE}.
##' @param Boundary.knots Boundary points at which to anchor the B-spline basis.
##' By default, they are the range of the non-\code{NA} data.  If both
##' \code{knots} and \code{Boundary.knots} are supplied, the basis parameters
##' do not depend on \code{x}. Data can extend beyond \code{Boundary.knots}.
##' @param ... Optional arguments for future usage.
##' @return A matrix of dimension \code{length(x)} by
##' \code{df = degree + length(knots)} (plus one if intercept is included).
##' Attributes that correspond to the arguments specified are returned
##' for usage for \code{\link{predict.bSpline2}}.
##' @examples
##' library(graphics)
##' x <- seq(0, 1, by = 0.01)
##' knots <- c(0.3, 0.5, 0.6)
##' bsMat <- bSpline(x, knots = knots, degree = 0, intercept = TRUE)
##' matplot(x, bsMat, type = "l", ylab = "B-spline basis with degree zero")
##' abline(v = knots, lty = 2, col = "gray")
##' @seealso
##' \code{\link{predict.bSpline2}} for evaluation at given (new) values;
##' \code{\link{ibs}} for integral of B-spline basis;
##' \code{\link{mSpline}} for M-spline basis;
##' \code{\link{iSpline}} for I-spline basis;
##' \code{\link{cSpline}} for C-spline basis.
##' @importFrom splines bs
##' @importFrom stats stepfun
##' @export
bSpline <- function(x, df = NULL, knots = NULL, degree = 3, intercept = FALSE,
                    Boundary.knots = range(x), ...) {

    ## check and reformat 'degree'
    if ((degree <- as.integer(degree)) < 0)
        stop("'degree' must be a nonnegative integer.")

    ## call splines::bs for non-zero degree
    if (degree > 0L) {
        out <- splines::bs(x = x, df = df, knots = knots, degree = degree,
                           intercept = intercept,
                           Boundary.knots = Boundary.knots)
        class(out) <- c("bSpline2", "basis", "matrix")
        return(out)
    }

    ## else degree is zero
    inputs <- pieceConst(x = x, df = df, knots = knots)
    knots <- inputs$knots
    ## potentially, df is a bad name since df is also a function in stats
    df <- inputs$df

    ## check whether any of x is outside of the boundary knots
    Boundary.knots <- sort(Boundary.knots)
    if (any(x < Boundary.knots[1] | x > Boundary.knots[2]))
        warning(paste("some 'x' values beyond boundary knots",
                      "may cause ill-conditioned bases"))

    ## piecewise constant basis
    augKnots <- c(Boundary.knots[1], knots, Boundary.knots[2])
    bsMat <- sapply(seq_len(df), function (i) {
        foo <- stats::stepfun(augKnots[i: (i + 1)], c(0L, 1L, 0L))
        foo(x)
    })

    ## include intercept or not
    if (! intercept)
        bsMat <- bsMat[, - 1L, drop = FALSE]

    ## add colnames for consistency with bs returns
    colnames(bsMat) <- as.character(seq_len(df  - as.integer(! intercept)))

    ## on attributes
    a <- list(degree = degree,
              knots = if (is.null(knots)) numeric(0L) else knots,
              Boundary.knots = Boundary.knots,
              intercept = intercept)
    attributes(bsMat) <- c(attributes(bsMat), a)
    class(bsMat) <- c("bSpline2", "basis", "matrix")
    bsMat
}


### internal function ==========================================================
##' @importFrom stats quantile
pieceConst <- function (x, df, knots) {
    ind <- (is.null(df) + 1) * is.null(knots) + 1
    ## ind == 1: knots is not NULL; df <- length(knots) + 1
    ## ind == 2: df is not NULL, while knots is NULL; number of piece <- df
    ## ind == 3: both df and knots are NULL; one-piece constant, df <- 1
    df <- switch(ind, length(knots) + 1L, as.integer(df), 1L)
    if (ind > 1) {
        tknots <- df + 1L
        quans <- seq.int(from = 0, to = 1,
                         length.out = tknots)[- c(1L, tknots)]
        knots <- as.numeric(stats::quantile(x, quans))
    }
    list(df = df, knots = knots)
}
