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
helper_update <- function(object, ..., fun, key_attr)
{
    dot_list <- list(...)
    if (length(dot_list) == 0L) {
        return(object)
    }
    check_attr(object, key_attr)
    if (is.null(dot_list$df)) {
        call_list <- modify_list(attributes(object), dot_list)
    } else {
        call_list <- modify_list(pred_attr(object, except = "knots"), dot_list)
    }
    do.call(fun, call_list)
}


##' @rdname update
##' @export
update.bSpline2 <- function(object, ...)
{
    helper_update(
        object,
        ...,
        fun = bSpline,
        key_attr = c("x", "degree", "knots", "Boundary.knots",
                     "intercept", "derivs", "integral")
    )
}


##' @rdname update
##' @export
update.dbs <- function(object, ...)
{
    helper_update(
        object,
        ...,
        fun = dbs,
        key_attr = c("x", "degree", "knots", "Boundary.knots",
                     "intercept", "derivs")
    )
}


##' @rdname update
##' @export
update.ibs <- function(object, ...)
{
    helper_update(
        object,
        ...,
        fun = ibs,
        key_attr = c("x", "degree", "knots", "Boundary.knots", "intercept")
    )
}


##' @rdname update
##' @export
update.mSpline <- function(object, ...)
{
    helper_update(
        object,
        ...,
        fun = mSpline,
        key_attr = c("x", "degree", "knots", "Boundary.knots",
                     "intercept", "derivs", "integral", "periodic")
    )
}


##' @rdname update
##' @export
update.iSpline <- function(object, ...)
{
    helper_update(
        object,
        ...,
        fun = iSpline,
        key_attr = c("x", "degree", "knots", "Boundary.knots",
                     "intercept", "derivs")
    )
}


##' @rdname update
##' @export
update.cSpline <- function(object, ...)
{
    helper_update(
        object,
        ...,
        fun = cSpline,
        key_attr = c("x", "degree", "knots", "Boundary.knots",
                     "intercept", "derivs", "scale")
    )
}


##' @rdname update
##' @export
update.bernsteinPoly <- function(object, ...)
{
    helper_update(
        object,
        ...,
        fun = bernsteinPoly,
        key_attr = c("x", "degree", "Boundary.knots",
                     "intercept", "derivs", "integral")
    )
}


##' @rdname update
##' @export
update.naturalSpline <- function(object, ...)
{
    helper_update(
        object,
        ...,
        fun = naturalSpline,
        key_attr = c("x", "knots", "Boundary.knots",
                     "intercept", "derivs", "integral")
    )
}
