################################################################################
##
##   R package splines2 by Wenjie Wang and Jun Yan
##   Copyright (C) 2016-2019
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


##' Print Out a Spline Basis Matrix
##'
##' \code{Print} methods that simply print out the spline basis matrix without
##' unnecessary attributes.
##'
##' @name print
##' @param x Objects of class \code{bSpline2}, \code{ibs}, \code{mSpline},
##' \code{iSpline}, or \code{cSpline}, etc.
##' @param ... Optional argument for future usage.
##'
##' @return Object input.
NULL


##' @rdname print
##' @export
print.bSpline2 <- function(x, ...) {
    print.default(tidyAttr(x, ...))
    invisible(x)
}


##' @rdname print
##' @export
print.ibs <- function(x, ...) {
    print.default(tidyAttr(x, ...))
    invisible(x)
}


##' @rdname print
##' @export
print.dbs <- function(x, ...) {
    print.default(tidyAttr(x, ...))
    invisible(x)
}


##' @rdname print
##' @export
print.mSpline <- function(x, ...) {
    print.default(tidyAttr(x, ...))
    invisible(x)
}

##' @rdname print
##' @export
print.iSpline <- function(x, ...) {
    print.default(tidyAttr(x, ...))
    invisible(x)
}


##' @rdname print
##' @export
print.cSpline <- function(x, ...) {
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
