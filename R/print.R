##
## R package splines2 by Wenjie Wang and Jun Yan
## Copyright (C) 2016-2021
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

##' @export
print.bSpline2 <- function(x, ...) {
    print.default(tidyAttr(x, ...))
    invisible(x)
}

##' @export
print.ibs <- function(x, ...) {
    print.default(tidyAttr(x, ...))
    invisible(x)
}

##' @export
print.dbs <- function(x, ...) {
    print.default(tidyAttr(x, ...))
    invisible(x)
}

##' @export
print.mSpline <- function(x, ...) {
    print.default(tidyAttr(x, ...))
    invisible(x)
}

##' @export
print.iSpline <- function(x, ...) {
    print.default(tidyAttr(x, ...))
    invisible(x)
}

##' @export
print.cSpline <- function(x, ...) {
    print.default(tidyAttr(x, ...))
    invisible(x)
}

##' @export
print.bernsteinPoly <- function(x, ...) {
    print.default(tidyAttr(x, ...))
    invisible(x)
}

##' @export
print.naturalSpline <- function(x, ...) {
    print.default(tidyAttr(x, ...))
    invisible(x)
}

### internal function ==========================================================
## remove all attributes but dim and dimnames
tidyAttr <- function(x, ...) {
    dimen <- attr(x, "dim")
    dimenName <- attr(x, "dimnames")
    attributes(x) <- NULL
    attr(x, "dim") <- dimen
    attr(x, "dimnames") <- dimenName
    x
}
