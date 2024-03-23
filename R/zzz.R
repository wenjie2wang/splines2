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

## set default options
splines2_default_options <- list(
    splines2.warn.outside = TRUE
)

## set options for splines2
.onLoad <- function(libname, pkgname)
{
  op <- options()
  toset <- ! names(splines2_default_options) %in% names(op)
  if (any(toset)) {
      options(splines2_default_options[toset])
  }
  invisible(NULL)
}
