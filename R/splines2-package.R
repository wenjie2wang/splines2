##
## R package splines2 by Wenjie Wang and Jun Yan
## Copyright (C) 2016-2025
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

##' splines2: Regression Spline Functions and Classes
##'
##' This package provides functions to construct basis matrices of
##' \itemize{
##' \item B-splines
##' \item M-splines
##' \item I-splines
##' \item convex splines (C-splines)
##' \item periodic splines
##' \item natural cubic splines
##' \item generalized Bernstein polynomials
##' \item along with their integrals (except C-splines) and derivatives
##'     of given order by closed-form recursive formulas
##' }
##'
##' In addition to the R interface, it also provides a C++ header-only library
##' integrated with \pkg{Rcpp}, which allows the construction of spline basis
##' functions directly in C++ with the help of \pkg{Rcpp} and
##' \pkg{RcppArmadillo}.  Thus, it can also be treated as one of the \pkg{Rcpp*}
##' packages.  A toy example package that uses the C++ interface is available at
##' <https://github.com/wenjie2wang/example-pkg-Rcpp-splines2>.
##'
##' The package \pkg{splines2} is intended to be a user-friendly supplement to
##' the base package \pkg{splines}.  The trailing number two in the package name
##' means "too" (and by no means refers to the generation two).  See Wang and
##' Yan (2021) for details and illustrations of how the package can be applied
##' to shape-restricted regression.
##'
##' @references
##'
##' Wang, W., & Yan, J. (2021). Shape-restricted regression splines with R
##' package \pkg{splines2}. \emph{Journal of Data Science}, 19(3), 498--517.
##'
##' @importFrom Rcpp sourceCpp
##' @useDynLib splines2
##'
##' @name splines2
##' @aliases splines2-package
##' @keywords internal
"_PACKAGE"
