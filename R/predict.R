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

##' Compute Spline Function for Given Coefficients
##'
##' Returns the spline function (with the specified coefficients) or evaluate
##' the basis functions at the specified \code{x} if the coefficients are not
##' specified.
##'
##' @name predict
##' @param object Spline objects produced by the \code{splines2} package.
##' @param newx The \code{x} values at which evaluations are required.  If it is
##'     \code{NULL} (by default), the original \code{x} used to create the
##'     spline object will be used.
##' @param coef A numeric vector specifying the coefficients of the spline basis
##'     functions.  If it is \code{NULL} (by default), the spline basis
##'     functions will be returned.  Otherwise, the resulting spline function
##'     will be returned.
##' @param ... Other options passed to the corresponding function that
##'     constructs the input \code{object}.  For example, the additional options
##'     will be passed to \code{bSpline()} for a \code{BSpline} object.
##'
##' @return The function returns the spline basis functions with the new values
##'     of \code{x} if \code{coef} is not specified.  Otherwise, the function
##'     returns the resulting spline function (or its derivative if
##'     \code{derivs} is specified as a positive integer through \code{...}).
##'
##' @example inst/examples/ex-predict.R
##'
##' @importFrom stats predict
NULL


## the helper function for predict
helper_predict <- function(object, newx = NULL, coef = NULL,
                           ..., .FUN, .KEY_ATTR)
{
    dot_list <- list(...)
    if (is.null(newx)) {
        if (is.null(coef) && length(dot_list) == 0) {
            return(object)
        }
        newx <- attr(object, "x")
    }
    ## checks if key attributes still exist
    check_attr(object, .KEY_ATTR)
    default_call <- c(list(x = newx), pred_attr(object))
    call_list <- modify_list(default_call, dot_list)
    out <- do.call(.FUN, call_list)
    if (is.null(coef)) {
        return(out)
    }
    ## else
    as.numeric(out %*% coef)
}


##' @rdname predict
##' @export
predict.BSpline <- function(object, newx = NULL, coef = NULL, ...)
{
    helper_predict(
        object = object,
        newx = newx,
        coef = coef,
        ...,
        .FUN = bSpline,
        .KEY_ATTR = c("x", "degree", "knots", "Boundary.knots",
                      "intercept", "periodic", "derivs", "integral")
    )
}

##' @rdname predict
##' @export
predict.MSpline <- function(object, newx = NULL, coef = NULL, ...)
{
    helper_predict(
        object = object,
        newx = newx,
        coef = coef,
        ...,
        .FUN = mSpline,
        .KEY_ATTR = c("x", "degree", "knots", "Boundary.knots",
                      "intercept", "periodic", "derivs", "integral")
    )
}


##' @rdname predict
##' @export
predict.ISpline <- function(object, newx = NULL, coef = NULL, ...)
{
    helper_predict(
        object = object,
        newx = newx,
        coef = coef,
        ...,
        .FUN = iSpline,
        .KEY_ATTR = c("x", "degree", "knots", "Boundary.knots",
                      "intercept", "derivs")
    )
}


##' @rdname predict
##' @export
predict.CSpline <- function(object, newx = NULL, coef = NULL, ...)
{
    helper_predict(
        object = object,
        newx = newx,
        coef = coef,
        ...,
        .FUN = cSpline,
        .KEY_ATTR = c("x", "degree", "knots", "Boundary.knots",
                      "intercept", "derivs", "scale")
    )
}


##' @rdname predict
##' @export
predict.BernsteinPoly <- function(object, newx = NULL, coef = NULL, ...)
{
    helper_predict(
        object = object,
        newx = newx,
        coef = coef,
        ...,
        .FUN = bernsteinPoly,
        .KEY_ATTR = c("x", "degree", "Boundary.knots",
                      "intercept", "derivs", "integral")
    )
}


##' @rdname predict
##' @export
predict.NaturalSpline <- function(object, newx = NULL, coef = NULL, ...)
{
    helper_predict(
        object = object,
        newx = newx,
        coef = coef,
        ...,
        .FUN = naturalSpline,
        .KEY_ATTR = c("x", "knots", "Boundary.knots",
                      "intercept", "derivs", "integral")
    )
}
