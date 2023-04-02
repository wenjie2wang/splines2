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

##' Visualize Spline Basis Functions
##'
##' Plot spline basis functions by lines in different colors.
##'
##' This function is intended to quickly visualize the spline basis functions.
##'
##' @param x A \code{splines2} object.
##' @param y An argument that is not used.
##' @param mark_knots A character vector specifying if knot placement should be
##'     indicated by vertical lines.
##' @param ... Additional arguments (other than \code{x} and \code{y}) that
##'     would be passed to \code{matplot()}.
##'
##' @importFrom graphics matplot abline
##' @export
plot.splines2 <- function(x, y,
                          mark_knots = c("none", "internal", "boundary", "all"),
                          ...)
{
    dots <- list(...)
    x_ <- attr(x, "x")
    default_args <- list(type = "l", xlab = "x", ylab = "")
    call_args <- modify_list(default_args, dots)
    call_args$x <- x_
    call_args$y <- x
    do.call(graphics::matplot, call_args)
    mark_knots <- match.arg(mark_knots,
                            choices = c("none", "internal", "boundary", "all"))
    ## prepare for marks
    if (mark_knots != "none" && isTRUE(attr(x, "periodic"))) {
        ik <- knots(x, "internal")
        bk <- knots(x, "boundary")
        dist_bk <- bk[2] - bk[1]
        range_x <- range(x_)
        k1 <- (bk[1] - range_x[1]) %/% dist_bk
        k2 <- (range_x[2] - bk[1]) %/% dist_bk
        shifts <- seq.int(k1, k2) * dist_bk
        ik <- unique(as.numeric(sapply(shifts, function(a) ik + a)))
        bk <- unique(as.numeric(sapply(shifts, function(a) bk + a)))
    }
    if (mark_knots == "internal" || mark_knots == "all") {
        graphics::abline(v = ik, col = "gray70", lty = 3)
    }
    if (mark_knots == "boundary" || mark_knots == "all") {
        graphics::abline(v = bk, col = "gray50", lty = 2)
    }
    invisible(x)
}
