##
## R package splines2 by Wenjie Wang and Jun Yan
## Copyright (C) 2016-2024
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

### some simple internal functions ===========================================

## wrap messages and keep proper line length
wrapMessages <- function(..., strwrap.args = list()) {
    x <- paste(...)
    wrap_x <- do.call(strwrap, c(list(x = x), strwrap.args))
    paste(wrap_x, collapse = "\n")
}

## convert null to numeric(0)
null2num0 <- function(x) {
    if (is.null(x)) {
        return(numeric(0))
    }
    x
}

## check key attributions
check_attr <- function(x, attrs = c("x", "degree", "knots",
                                    "Boundary.knots", "intercept"))
{
    idx <- ! attrs %in% names(attributes(x))
    if (any(idx)) {
        stop(wrapMessages(
            sprintf("Missing attributes: %s.",
                    paste(attrs[idx], collapse = ", "))
        ), call. = FALSE)
    }
    invisible()
}

## return most of the attributes except x and class
pred_attr <- function(x, except = c("x", "class", "dimnames")) {
    out <- attributes(x)
    out[! names(out) %in% except]
}

## simplified version of utils::modifyList with keep.null = TRUE
modify_list <- function (x, val) {
    stopifnot(is.list(x), is.list(val))
    xnames <- names(x)
    vnames <- names(val)
    vnames <- vnames[nzchar(vnames)]
    for (v in vnames) {
        x[v] <- if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]]))
                    list(modify_list(x[[v]], val[[v]]))
                else val[v]
    }
    x
}
