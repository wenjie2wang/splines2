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
helper_makepredictcall <- function(var, call, .FUN, .KEY_ATTR)
{
    fun_symbol <- substitute(.FUN)
    fun_name <- as.character(fun_symbol)

    ## the following checking seems to be unnecessary?
    ## if (as.character(call)[1L] == fun_name ||
    ##     (is.call(call) && identical(eval(call[[1L]]), eval(fun_symbol)))) {
    ##     at <- attributes(var)[key_attr]
    ##     call <- call[1L:2L]
    ##     call[names(at)] <- at
    ## }

    ## throw warnings instead
    res <- tryCatch(check_attr(var, .KEY_ATTR), error = function(e) e)
    if (inherits(res, "error")) {
        warning(res, call. = FALSE)
    } else {
        at <- attributes(var)[.KEY_ATTR]
        call <- call[1L:2L]
        call[names(at)] <- at
    }
    ## return
    call
}

##' @export
makepredictcall.BSpline <- function(var, call)
{
    helper_makepredictcall(
        var, call,
        .FUN = bSpline,
        .KEY_ATTR = c("degree", "knots", "Boundary.knots", "intercept",
                      "periodic", "derivs", "integral")
    )
}

##' @export
makepredictcall.NaturalSpline <- function(var, call)
{
    helper_makepredictcall(
        var, call,
        .FUN = naturalSpline,
        .KEY_ATTR = c("knots", "Boundary.knots", "trim",
                      "intercept", "derivs", "integral")
    )
}

##' @export
makepredictcall.NaturalSplineK <- function(var, call)
{
    helper_makepredictcall(
        var, call,
        .FUN = nsk,
        .KEY_ATTR = c("knots", "Boundary.knots", "trim",
                      "intercept", "derivs", "integral")
    )
}

##' @export
makepredictcall.MSpline <- function(var, call)
{
    helper_makepredictcall(
        var, call,
        .FUN = mSpline,
        .KEY_ATTR = c("degree", "knots", "Boundary.knots", "intercept",
                      "periodic", "derivs", "integral")
    )
}

##' @export
makepredictcall.ISpline <- function(var, call)
{
    helper_makepredictcall(
        var, call,
        .FUN = iSpline,
        .KEY_ATTR = c("degree", "knots", "Boundary.knots", "intercept",
                      "derivs")
    )
}

##' @export
makepredictcall.CSpline <- function(var, call)
{
    helper_makepredictcall(
        var, call,
        .FUN = cSpline,
        .KEY_ATTR = c("degree", "knots", "Boundary.knots", "intercept",
                      "derivs", "scale")
    )
}

##' @export
makepredictcall.BernsteinPoly <- function(var, call)
{
    helper_makepredictcall(
        var, call,
        .FUN = bernsteinPoly,
        .KEY_ATTR = c("degree", "Boundary.knots", "intercept",
                      "derivs", "integral")
    )
}
