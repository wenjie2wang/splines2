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


#' Bernstein Polynomial Basis
#'
#' Generate the Bernstein spline basis matrix given a degree of freedom
#'
#' @param x The predictor variable. Missing values are allowed.
#' @param df Degrees of freedom.
#' ## @param Boundary.knots Boundary points at wich to anchor the Bernstein basis.
#' @param ... Optional arguments for future usage.
#' @return A matrix of dimension \code{length(x)} by \code{df}.
#' TODO: Attributes of derivatives, integral, etc.
#' @references
#' https://people.sc.fsu.edu/~jburkardt/c_src/bernstein_polynomial/bernstein_polynomial.html
#'
#' @examples
#'
#'
bernsteinSpline <- function(x, df, ...) {
    ## z <- (x - Boundary.knots[1]) / diff(Boundary.knots)
    stopifnot(all(x >= 0 & x <= 1))
    z <- x
    n <- df - 1
    idx <- 0:n
    zmat <- sapply(idx, function(i) choose(n, i) *  z^i * (1 - z)^(n - i) * df)
    xmat <- zmat # * diff(Boundary.knots) + Boundary.knots[1]
    ## TODO: set attributes for derivatives
    class(xmat) <- c("bernsteinSpline", "basis", "matrix")
    xmat
}

## Derivatives
evalA <- function(z, Btahat) {
  df <- length(Btahat)
  xbasis <- bernsteinSpline(z, df)
  c(xbasis %*% Btahat)
}


d1A <- function(z, Btahat) {
  df <- length(Btahat)
  xbasis <- bernsteinSpline(z, df - 1)
  dBtahat <- diff(Btahat)
  N <- df - 1
  N * c(xbasis %*% dBtahat)
}


d2A <- function(z, Btahat) {
  df <- length(Btahat)
  xbasis <- bernsteinSpline(z, df - 2)
  ddBtahat <- diff(diff(Btahat))
  N <- df - 1
  N * (N - 1) * c(xbasis %*% ddBtahat)
}
