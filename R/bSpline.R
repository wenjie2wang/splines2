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


##' B-Spline Basis for Polynomial Splines
##'
##' This function generates the B-spline basis matrix for a polynomial spline.
##'
##' It is an augmented function of \code{\link[splines]{bs}} in package
##' \code{splines} for B-spline basis that allows piecewise constant (close
##' on the left, open on the right) spline basis with zero degree. When the
##' argument \code{degree} is greater than zero, it internally calls
##' \code{\link[splines]{bs}} and generates a basis matrix for representing
##' the family of piecewise polynomials with the specified interior knots and
##' degree, evaluated at the values of \code{x}.  The function has the same
##' arguments with \code{\link[splines]{bs}} for ease usage.
##'
##' @usage
##' bSpline(x, df = NULL, knots = NULL, degree = 3L, intercept = FALSE,
##'         Boundary.knots = range(x, na.rm = TRUE), ...)
##'
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
##'
##' @return A matrix of dimension \code{length(x)} by
##' \code{df = degree + length(knots)} (plus one if intercept is included).
##' Attributes that correspond to the arguments specified are returned
##' for usage of other functions in this package.
##' @examples
##' library(splines2)
##' x <- seq(0, 1, by = 0.01)
##' knots <- c(0.3, 0.5, 0.6)
##' bsMat <- bSpline(x, knots = knots, degree = 0, intercept = TRUE)
##'
##' library(graphics)
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
bSpline <- function(x, df = NULL, knots = NULL, degree = 3L, intercept = FALSE,
                    Boundary.knots = range(x, na.rm = TRUE), ...)
{
    ## check and reformat 'degree'
    if ((degree <- as.integer(degree)) < 0)
        stop("'degree' must be a nonnegative integer.")

    ## sort and remove possible NA's in internal knots if exist
    if (length(knots)) {
        knots <- sort(knots)
        tmp <- is.na(knots)
        if (any(tmp)) {
            omit <- seq_along(knots)[tmp]
            knots <- knots[- omit]
        }
    }

    ## take care of possible NA's in `x` for `Boundary.knots`
    nax <- is.na(x)
    if (all(nax))
        stop("'x' cannot be all NA's!")
    ## if (missing(Boundary.knots) || is.null(Boundary.knots) ||
    ##     any(is.na(Boundary.knots)))
    ##     Boundary.knots <- range(x[! nax])

    ## call splines::bs for non-zero degree
    if (degree > 0L) {
        out <- splines::bs(x = x, df = df, knots = knots, degree = degree,
                           intercept = intercept,
                           Boundary.knots = Boundary.knots)
        attr(out, "x") <- x
        class(out) <- c("matrix", "bSpline2")
        return(out)
    }

    ## else degree is zero
    xx <- x
    ## remove NA's in x
    if (nas <- any(nax))
        xx <- x[! nax]

    ## prepare inputs for piecewise constant bases
    inputs <- pieceConst(x = xx, df = df, knots = knots)
    knots <- inputs$knots
    ## potentially, df is a bad name since df is also a function in stats
    df <- inputs$df

    ## check whether any of x is outside of the boundary knots
    Boundary.knots <- sort(Boundary.knots)
    if (any(xx < Boundary.knots[1L] | xx > Boundary.knots[2L]))
        warning(paste("Some 'x' values beyond boundary knots",
                      "may cause ill-conditioned bases."))

    ## piecewise constant basis
    augKnots <- c(Boundary.knots[1L], knots, Boundary.knots[2L])
    bsMat <- sapply(seq_len(df), function (i) {
        foo <- stats::stepfun(augKnots[i: (i + 1L)], c(0L, 1L, 0L))
        foo(xx)
    })

    ## make sure bsMat is a matrix
    if (! is.matrix(bsMat))
        bsMat <- matrix(bsMat, nrow = length(xx))

    ## include intercept or not
    if (! intercept) {
        if (length(knots))
            bsMat <- bsMat[, - 1L, drop = FALSE]
        else
            stop(paste("'intercept' has to be 'TRUE'",
                       "for one-piece const basis."))
    }

    ## keep NA's as is
    if (nas) {
        nmat <- matrix(NA, length(nax), ncol(bsMat))
        nmat[! nax, ] <- bsMat
        bsMat <- nmat
    }

    ## add colnames for consistency with bs returns
    colnames(bsMat) <- as.character(seq_len(df - as.integer(! intercept)))

    ## on attributes
    tmp <- list(degree = degree,
                knots = if (is.null(knots)) numeric(0L) else knots,
                Boundary.knots = Boundary.knots,
                intercept = intercept, x = x)
    attributes(bsMat) <- c(attributes(bsMat), tmp)
    class(bsMat) <- c("matrix", "bSpline2")
    bsMat
}


### internal function ==========================================================
##' @importFrom stats quantile
pieceConst <- function (x, df, knots)
{
    if (any(is.na(knots)))
        stop("'knots' cannot contain any NA.")
    if (! is.null(knots))
        knots <- sort(knots)
    ind <- (is.null(df) + 1L) * is.null(knots) + 1L
    ## ind == 1: knots is not NULL; df <- length(knots) + 1
    ## ind == 2: df is not NULL, while knots is NULL; number of piece <- df
    ## ind == 3: both df and knots are NULL; one-piece constant, df <- 1
    df <- switch(ind, length(knots) + 1L, as.integer(df), 1L)
    if (df < 1L) {
        df <- 1L
        warning("'df' specified was too small; used 1 instead.")
    }
    if (ind > 1) {
        tknots <- df + 1L
        quans <- seq.int(from = 0, to = 1,
                         length.out = tknots)[- c(1L, tknots)]
        knots <- as.numeric(stats::quantile(x, quans))
    }
    list(df = df, knots = knots)
}
