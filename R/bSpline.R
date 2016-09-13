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


##'
##' @importFrom splines bs
##' @importFrom stats stepfun quantile
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
        return(out)
    }

    ## else degree is zero
    inputs <- pieceConst(x = x, df = df, knots = knots)
    knots <- inputs$knots
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
    if (! intercept)
        bsMat <- bsMat[, - 1L, drop = FALSE]
    a <- list(degree = degree,
             knots = if (is.null(knots)) numeric(0L) else knots,
             Boundary.knots = Boundary.knots,
             intercept = intercept)
    attributes(bsMat) <- c(attributes(bsMat), a)
    class(bsMat) <- c("bs", "basis", "matrix")
    bsMat
}


### internal function ==========================================================
pieceConst <- function (x, df, knots) {
    ind <- (is.null(df) + 1) * is.null(knots) + 1
    ## ind == 1: knots is not NULL; df <- length(knots) + 1
    ## ind == 2: df is not NULL, while knots is NULL; number of piece <- df
    ## ind == 3: both df and knots are NULL; one-piece constant, df <- 1
    df <- switch(ind, length(knots) + 1L, as.integer(df), 1L)
    if (ind > 1) {
        tknots <- df + 1L
        quans <- seq.int(from = 0, to = 1,
                        length.out = tknots)[-c(1L, tknots)]
        knots <- as.numeric(stats::quantile(x, quans))
    }
    ## return
    list(df = df, knots = knots)
}
