##
## R package splines2 by Wenjie Wang and Jun Yan
## Copyright (C) 2016-2022
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

##' Extract Knots from the Given Object
##'
##' Methods for the generic function \code{knots} from the \pkg{stats} package
##' to obtain internal or boundary knots from the objects produced by this
##' package.
##'
##' @name knots
##'
##' @param Fn An \code{splines2} object produced by this package.
##' @param type A character vector of length one indicating the type of knots to
##'     return.  The available choices are \code{"internal"} for internal knots
##'     and \code{"Boundary"} for boundary knots.
##' @param ... Optional arguments that are not used now.
##'
##' @return A numerical vector.
##'
##' @example inst/examples/ex-knots.R
##'
##' @importFrom stats knots
NULL


## the default method
knots_ <- function(Fn, type = c("internal", "boundary"), ...) {
    type <- match.arg(type, choices = c("internal", "boundary"))
    if (type == "internal") {
        attr(Fn, "knots")
    } else {
        attr(Fn, "Boundary.knots")
    }
}


##' @rdname knots
##' @export
knots.bSpline2 <- knots_

##' @rdname knots
##' @export
knots.dbs <- knots_

##' @rdname knots
##' @export
knots.ibs <- knots_

##' @rdname knots
##' @export
knots.mSpline <- knots_

##' @rdname knots
##' @export
knots.iSpline <- knots_

##' @rdname knots
##' @export
knots.cSpline <- knots_

##' @rdname knots
##' @export
knots.bernsteinPoly <- knots_

##' @rdname knots
##' @export
knots.naturalSpline <- knots_
