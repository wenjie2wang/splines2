[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/splines2)](http://cran.r-project.org/package=splines2)

# splines2

The R package **splines2** is a complementary package on splines providing
functions constructing B-spline, monotone spline (M-spline) and its integral
(I-spline), convex spline (C-spline), and integral of B-spline basis. Piecewise
constant basis is allowed for B-spline and M-spline basis.


## Development

The latest version of package is under development at
[GitHub](https://github.com/wenjie2wang/splines2) in branch 'dev'.


## Installation of Stable Version

You can install the stable version from [CRAN](https://cran.r-project.org/package=splines2).

```r
install.packages("splines2", dependencies = TRUE)
```


## Basic Usage

```r
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

The R package splines2 is free software: You can redistribute it and/or
modify it under the terms of the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License, or
any later version (at your option).
See the [GNU General Public License](http://www.gnu.org/licenses/) for details.

The R package splines2 is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
