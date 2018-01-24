################################################################################
##
##   R package splines2 by Wenjie Wang and Jun Yan
##   Copyright (C) 2016-2018
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


### some trivial internal functions ============================================
## wrap messages and keep proper line length
wrapMessages <- function(..., strwrap.args = list()) {
    x <- paste(...)
    wrap_x <- do.call(strwrap, c(list(x = x), strwrap.args))
    paste(wrap_x, collapse = "\n")
}

## warning if x contains NA (or NaN)
na_warning <- function(x, sub_env = c("current", "parent", "grandparent"),
                       num_grandparent = 2L, ...)
{
    sub_env <- switch(
        match.arg(sub_env),
        "current" = environment(),
        "parent" = parent.frame(),
        "grandparent" = parent.frame(num_grandparent)
    )
    objName = deparse(substitute(x, sub_env))
    if (anyNA(x))
        warning(wrapMessages(
            sprintf("Found `NA` values in `%s`.", objName)
        ), call. = FALSE)
    invisible(x)
}

## stop if x contains NA (or NaN)
na_stop <- function(x, sub_env = c("current", "parent", "grandparent"),
                       num_grandparent = 2L, ...)
{
    sub_env <- switch(
        match.arg(sub_env),
        "current" = environment(),
        "parent" = parent.frame(),
        "grandparent" = parent.frame(num_grandparent)
    )
    objName = deparse(substitute(x, sub_env))
    if (anyNA(x))
        stop(wrapMessages(
            sprintf("Found `NA` values in `%s`.", objName)
        ), call. = FALSE)
    invisible(x)
}

## is x a numeric matrix (optionally of nRow rows and nCol columns)
isNumMatrix <- function(x, nRow = NULL, nCol = NULL,
                        warn_na = TRUE, error_na = ! warn_na,
                        sub_env = "parent", ...)
{
    out <- is.numeric(x) && is.matrix(x)
    if (out) {
        nDim <- dim(x)
        if (! is.null(nRow)) out <- out && nDim[1L] == nRow
        if (! is.null(nCol)) out <- out && nDim[2L] == nCol
        if (error_na) na_stop(x, sub_env = sub_env, ...)
        if (warn_na) na_warning(x, sub_env = sub_env, ...)
    }
    out
}
