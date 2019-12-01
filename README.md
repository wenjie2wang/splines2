# splines2

[![CRAN_Status_Badge][r-pkg-badge]][cran-url]
[![Build Status][travis-master]][travis]
[![AppVeyor Build Status][appveyor-master]][appveyor]
[![codecov][codecov-master]][codecov]
[![Downloads from the RStudio CRAN mirror][cranlog-badge]][cran-url]

The R package **splines2** is a supplementary package on splines providing
functions constructing B-splines, integral of B-splines, monotone splines
(M-splines) and its integral (I-splines), convex splines (C-splines), and their
derivatives of given order. Piecewise constant basis is allowed for B-spline and
M-spline basis.


## Installation of CRAN Version

You can install the released version from [CRAN][cran-url].

```R
install.packages("splines2")
```

## Getting Started

[Online documentation][pkgdown-url] provides reference for all functions and
contains package vignettes for
- A quick [demonstration][pkg-demo] for the basic usage of the main functions


## Development


The latest version of package is under development at [GitHub][github-url].  If
it is able to pass the building check by Travis CI, you may consider installing
it with the help of **remotes** by

```R
if (! require(remotes)) install.packages("remotes")
remotes::install_github("wenjie2wang/splines2")
```


## License

The R package **splines2** is free software: You can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or any later version
(at your option).  See the [GNU General Public License][gpl] for details.

The R package **splines2** is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.


[r-pkg-badge]: https://www.r-pkg.org/badges/version/splines2
[cranlog-badge]: https://cranlogs.r-pkg.org/badges/splines2
[cran-url]: https://CRAN.R-project.org/package=splines2
[travis]: https://travis-ci.org/wenjie2wang/splines2
[travis-master]: https://travis-ci.org/wenjie2wang/splines2.svg?branch=master
[github-url]: https://github.com/wenjie2wang/splines2
[pkgdown-url]: https://wenjie-stat.me/splines2
[pkg-demo]: https://wenjie-stat.me/splines2/articles/splines2-intro
[codecov]: https://codecov.io/gh/wenjie2wang/splines2
[codecov-master]: https://codecov.io/gh/wenjie2wang/splines2/branch/master/graph/badge.svg
[appveyor]: https://ci.appveyor.com/project/wenjie2wang/splines2
[appveyor-master]: https://ci.appveyor.com/api/projects/status/bvoso7nxchg1incb/branch/master?svg=true
[gpl]: https://www.gnu.org/licenses/
