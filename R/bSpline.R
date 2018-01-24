################################################################################
##
##   R package splines2 by Wenjie Wang and Jun Yan
##   Copyright (C) 2016-2018
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
##' \code{splines} for B-spline basis that allows piecewise constant (close on
##' the left, open on the right) spline basis with zero degree. When the
##' argument \code{degree} is greater than zero, it internally calls
##' \code{\link[splines]{bs}} and generates a basis matrix for representing the
##' family of piecewise polynomials with the specified interior knots and
##' degree, evaluated at the values of \code{x}.  The function has the same
##' arguments with \code{\link[splines]{bs}} for ease usage.
##'
##' @usage
##' bSpline(x, df = NULL, knots = NULL, degree = 3L, intercept = FALSE,
##'         Boundary.knots = range(x, na.rm = TRUE), ...)
##'
##' @param x The predictor variable.  Missing values are allowed and will be
##'     returned as they were.
##' @param df Degrees of freedom.  One can specify \code{df} rather than
##'     \code{knots}, then the function chooses "df - degree" (minus one if
##'     there is an intercept) knots at suitable quantiles of \code{x} (which
##'     will ignore missing values).  The default, \code{NULL}, corresponds to
##'     no inner knots, i.e., "degree - intercept". If \code{knots} was
##'     specified, \code{df} specified will be ignored.
##' @param knots The internal breakpoints that define the spline.  The default
##'     is \code{NULL}, which results in a basis for ordinary polynomial
##'     regression.  Typical values are the mean or median for one knot,
##'     quantiles for more knots.  See also \code{Boundary.knots}.
##' @param degree Non-negative integer degree of the piecewise polynomial. The
##'     default value is 3 for cubic splines. Zero degree is allowed for this
##'     function, which is the only difference compared with
##'     \code{\link[splines]{bs}} in package \code{splines}.
##' @param intercept If \code{TRUE}, an intercept is included in the basis;
##'     Default is \code{FALSE}.
##' @param Boundary.knots Boundary points at which to anchor the B-spline basis.
##'     By default, they are the range of the non-\code{NA} data.  If both
##'     \code{knots} and \code{Boundary.knots} are supplied, the basis
##'     parameters do not depend on \code{x}. Data can extend beyond
##'     \code{Boundary.knots}.
##' @param ... Optional arguments for future usage.
##'
##' @return A matrix of dimension \code{length(x)} by
##' \code{df = degree + length(knots)} (plus one if intercept is included).
##' Attributes that correspond to the arguments specified are returned
##' for usage of other functions in this package.
##' @examples
##' library(splines2)
##' x <- seq.int(0, 1, 0.01)
##' knots <- c(0.3, 0.5, 0.6)
##' bsMat <- bSpline(x, knots = knots, degree = 0, intercept = TRUE)
##'
##' library(graphics)
##' matplot(x, bsMat, type = "l", ylab = "Piecewise constant B-spline bases")
##' abline(v = knots, lty = 2, col = "gray")
##' @seealso
##' \code{\link{predict.bSpline2}} for evaluation at given (new) values;
##' \code{\link{dbs}}, \code{\link{deriv.bSpline2}} for derivatives;
##' \code{\link{ibs}} for integral of B-splines;
##' \code{\link{mSpline}} for M-splines;
##' \code{\link{iSpline}} for I-splines;
##' \code{\link{cSpline}} for C-splines.
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
    if (length(knots))
        knots <- sort.int(knots)

    ## take care of possible NA's in `x`
    nax <- is.na(x)
    if (all(nax))
        stop("The 'x' cannot be all NA's!")

    ## call splines::bs for non-zero degree
    if (degree) {
        out <- splines::bs(x = x, df = df, knots = knots,
                           degree = degree, intercept = intercept,
                           Boundary.knots = Boundary.knots)
        ## add "x" to attributes
        attr(out, "x") <- x
        ## throw out warning if any internal knot outside boundary.knots
        knots <- attr(out, "knots")
        Boundary.knots <- attr(out, "Boundary.knots")
        ## any internal knots placed outside of boundary knots?
        outside_knots <- (knots <= Boundary.knots[1L]) |
            (knots >= Boundary.knots[2L])
        if (any(outside_knots))
            warning(wrapMessages(
                "Some internal knots were not placed",
                "inside of boundary knots,",
                "which may cause \nill-conditioned bases!"
            ))
        ## update classes
        class(out) <- c("matrix", "bSpline2")
        return(out)
    }

    ## else degree is zero
    ## remove NA's in x
    xx <- if (nas <- any(nax)) x[! nax] else x

    ## check whether any of x is outside of the boundary knots
    outside_x <- rep(FALSE, length(xx))
    if (! missing(Boundary.knots)) {
        if (! is.numeric(Boundary.knots) || anyNA(Boundary.knots))
            stop(wrapMessages(
                "The 'Boundary.knots' has to be",
                "numeric vector of length 2",
                "with no missing value."
            ))
        if (length(Boundary.knots) > 2) {
            warning(wrapMessages(
                "Only the first two values",
                "in the 'Boundary.knots' were used."
            ))
            Boundary.knots <- Boundary.knots[seq_len(2L)]
        }
        Boundary.knots <- sort.int(Boundary.knots)
        outside_x <- (xx < Boundary.knots[1L]) | (xx > Boundary.knots[2L])
    }
    if (any(outside_x))
        warning(wrapMessages(
            "Some 'x' values beyond boundary knots",
            "may cause ill-conditioned bases!"
        ))

    ## prepare inputs for piecewise constant bases
    inputs <- pieceConst(x = xx[! outside_x], df = df, knots = knots,
                         Boundary.knots = Boundary.knots)
    knots <- inputs$knots
    ## potentially, df is a bad name since df is also a function in stats
    df <- inputs$df

    ## piecewise constant basis
    augKnots <- c(Boundary.knots[1L], knots, Boundary.knots[2L])
    bsMat <- sapply(seq_len(df), function (i) {
        foo <- stats::stepfun(augKnots[i: (i + 1L)], c(0L, 1L, 0L))
        foo(xx)
    })

    ## close on the right boundary knot for the last constant piece?
    ## if (any(rightX <- xx == Boundary.knots[2L]))
    ##     bsMat[rightX, df] <- 1

    ## make sure bsMat is a matrix
    if (! is.matrix(bsMat))
        bsMat <- matrix(bsMat, nrow = length(xx))

    ## include intercept or not
    if (! intercept) {
        if (length(knots))
            bsMat <- bsMat[, - 1L, drop = FALSE]
        else
            stop(wrapMessages(
                "The 'intercept' has to be 'TRUE'",
                "for one-piece const basis."
            ))
    }

    ## keep NA's as is
    if (nas) {
        nmat <- matrix(NA, length(nax), ncol(bsMat))
        nmat[! nax, ] <- bsMat
        bsMat <- nmat
    }

    ## add dimnames for consistency with bs returns
    row.names(bsMat) <- names(x)
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
pieceConst <- function (x, df, knots, Boundary.knots)
{
    ind <- (is.null(df) + 1L) * is.null(knots) + 1L
    ## ind == 1: knots is not NULL; df <- length(knots) + 1
    ## ind == 2: df is not NULL, while knots is NULL; number of piece <- df
    ## ind == 3: both df and knots are NULL; one-piece constant, df <- 1
    df0 <- switch(ind, length(knots) + 1L, as.integer(df), 1L)
    if (ind > 1L) {
        tknots <- df0 + 1L
        quans <- seq.int(from = 0, to = 1,
                         length.out = tknots)[- c(1L, tknots)]
        knots <- as.numeric(stats::quantile(x, quans))
    } else {
        ## any internal knots placed outside of boundary knots?
        outside_knots <- (knots <= Boundary.knots[1L]) |
            (knots >= Boundary.knots[2L])
        ## remove internal knots placed outside of boundary knots
        if (any(outside_knots)) {
            knots <- knots[! outside_knots]
            df0 <- df0 - sum(outside_knots)
            warning(wrapMessages(
                "Only internal knots placed inside",
                "boundary knots were considered."
            ), call. = FALSE)
        }
        if (! is.null(df) && df != df0)
            warning(wrapMessages(
                "The 'df' specified was not appropriate.",
                sprintf("Used 'df = %d' instead.", df0)
            ),  call. = FALSE)
    }
    list(df = df0, knots = knots)
}
