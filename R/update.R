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

##' Update Spline Basis Functions
##'
##' Update the knot placement, polynomial degree, and any other options
##' available when constructing the given spline object.
##'
##' @name update
##'
##' @inheritParams predict
##'
##' @param ... Other arguments passed to the corresponing constructor function.
##'
##' @return An updated object of the same class as the input object with the
##'     specified updates.
##'
##' @example inst/examples/ex-update.R
##'
##' @importFrom stats update
NULL


## the helper function for update
helper_update <- function(object, ..., .FUN, .KEY_ATTR)
{
    dot_list <- list(...)
    if (length(dot_list) == 0L) {
        return(object)
    }
    check_attr(object, .KEY_ATTR)
    exclude_attrs <- c("class", "dimnames")
    if (! is.null(dot_list$df)) {
        exclude_attrs <- c(exclude_attrs, "knots")
    }
    if (! is.null(dot_list$trim)) {
        exclude_attrs <- c(exclude_attrs, "Boundary.knots")
    }
    call_attr <- pred_attr(object, except = exclude_attrs)
    call_list <- modify_list(call_attr, dot_list)
    do.call(.FUN, call_list)
}


##' @rdname update
##' @export
update.BSpline <- function(object, ...)
{
    helper_update(
        object,
        ...,
        .FUN = bSpline,
        .KEY_ATTR = c("x", "degree", "knots", "Boundary.knots",
                      "intercept", "derivs", "integral", "periodic")
    )
}


##' @rdname update
##' @export
update.MSpline <- function(object, ...)
{
    helper_update(
        object,
        ...,
        .FUN = mSpline,
        .KEY_ATTR = c("x", "degree", "knots", "Boundary.knots",
                      "intercept", "derivs", "integral", "periodic")
    )
}


##' @rdname update
##' @export
update.ISpline <- function(object, ...)
{
    helper_update(
        object,
        ...,
        .FUN = iSpline,
        .KEY_ATTR = c("x", "degree", "knots", "Boundary.knots",
                      "intercept", "derivs")
    )
}


##' @rdname update
##' @export
update.CSpline <- function(object, ...)
{
    helper_update(
        object,
        ...,
        .FUN = cSpline,
        .KEY_ATTR = c("x", "degree", "knots", "Boundary.knots",
                      "intercept", "derivs", "scale")
    )
}


##' @rdname update
##' @export
update.BernsteinPoly <- function(object, ...)
{
    helper_update(
        object,
        ...,
        .FUN = bernsteinPoly,
        .KEY_ATTR = c("x", "degree", "Boundary.knots",
                      "intercept", "derivs", "integral")
    )
}


##' @rdname update
##' @export
update.NaturalSpline <- function(object, ...)
{
    helper_update(
        object,
        ...,
        .FUN = naturalSpline,
        .KEY_ATTR = c("x", "knots", "Boundary.knots", "trim",
                      "intercept", "derivs", "integral")
    )
}

##' @rdname update
##' @export
update.NaturalSplineK <- function(object, ...)
{
    helper_update(
        object,
        ...,
        .FUN = nsk,
        .KEY_ATTR = c("x", "knots", "Boundary.knots", "trim",
                      "intercept", "derivs", "integral")
    )
}
