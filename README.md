[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/splines2)](http://cran.r-project.org/package=splines2)
[![Build Status](https://travis-ci.org/wenjie2wang/splines2.svg?branch=master)](https://travis-ci.org/wenjie2wang/splines2)


# splines2

The R package **splines2** is a complementary package on splines providing
functions constructing B-spline, monotone spline (M-spline) and its integral
(I-spline), convex spline (C-spline), and integral of B-spline basis. Piecewise
constant basis is allowed for B-spline and M-spline basis.


## Installation of CRAN Version

You can install the stable version from [CRAN](https://cran.r-project.org/package=splines2).

```R
install.packages("splines2", dependencies = TRUE)
```


## Development

The latest version of package is under development at
[GitHub](https://github.com/wenjie2wang/splines2) in branch 'dev'.

[![Build Status](https://travis-ci.org/wenjie2wang/splines2.svg?branch=dev)](https://travis-ci.org/wenjie2wang/splines2)

You may install the latest version under development with the help of *devtools*
if it passed the building check by Travis CI.

```R
if (! require(devtools)) install.packages("devtools", dependencies = TRUE)
devtools::install_git("git://github.com/wenjie2wang/splines2.git", branch = "dev")
```


## Basic Usage

```R
help(pacakge = "splines2")
library(splines2)
## for integral of B-spline basis
?ibs
## for M-spline basis
?mSpline
## for I-spline basis
?iSpline
```

[package help manual](https://cran.r-project.org/web/packages/splines2/splines2.pdf)
is also available for details and demonstration.


## License

The R package splines2 is free software: You can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or any later version (at
your option).  See
the [GNU General Public License](http://www.gnu.org/licenses/) for details.

The R package splines2 is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.
