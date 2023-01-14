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

## reference: splines:::makepredictcall

##' @importFrom stats makepredictcall
helper_makepredictcall <- function(var, call, fun, key_attr)
{
    fun_symbol <- substitute(fun)
    fun_name <- as.character(fun_symbol)
    if (as.character(call)[1L] == fun_name ||
        (is.call(call) && identical(eval(call[[1L]]), eval(fun_symbol)))) {
        at <- attributes(var)[key_attr]
        call <- call[1L:2L]
        call[names(at)] <- at
    }
    call
}

##' @export
makepredictcall.bSpline2 <- function(var, call)
{
    helper_makepredictcall(
        var, call, fun = bSpline,
        key_attr = c("degree", "knots", "Boundary.knots", "intercept",
                     "derivs", "integral")
    )
}

##' @export
makepredictcall.naturalSpline <- function(var, call)
{
    helper_makepredictcall(
        var, call, fun = naturalSpline,
        key_attr = c("knots", "Boundary.knots", "intercept",
                     "derivs", "integral")
    )
}

##' @export
makepredictcall.dbs <- function(var, call)
{
    helper_makepredictcall(
        var, call, fun = dbs,
        key_attr = c("degree", "knots", "Boundary.knots", "intercept", "derivs")
    )
}

##' @export
makepredictcall.ibs <- function(var, call)
{
    helper_makepredictcall(
        var, call, fun = ibs,
        key_attr = c("degree", "knots", "Boundary.knots", "intercept")
    )
}

##' @export
makepredictcall.mSpline <- function(var, call)
{
    helper_makepredictcall(
        var, call, fun = mSpline,
        key_attr = c("degree", "knots", "Boundary.knots", "intercept",
                     "periodic", "derivs", "integral")
    )
}

##' @export
makepredictcall.iSpline <- function(var, call)
{
    helper_makepredictcall(
        var, call, fun = iSpline,
        key_attr = c("degree", "knots", "Boundary.knots", "intercept", "derivs")
    )
}

##' @export
makepredictcall.cSpline <- function(var, call)
{
    helper_makepredictcall(
        var, call, fun = cSpline,
        key_attr = c("degree", "knots", "Boundary.knots", "intercept",
                     "derivs", "scale")
    )
}

##' @export
makepredictcall.bernsteinPoly <- function(var, call)
{
    helper_makepredictcall(
        var, call, fun = bernsteinPoly,
        key_attr = c("degree", "Boundary.knots", "intercept",
                     "derivs", "integral")
    )
}
