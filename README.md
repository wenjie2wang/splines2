# splines2

The R package **splines2** is a complementary package on splines providing
functions constructing B-splines, monotone splines (M-splines) and its integral
(I-splines), convex splines (C-splines), and integral of B-splines. Piecewise
constant basis of degree zero is allowed for B-spline and M-spline basis.


## Installation of CRAN Version

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/splines2)](https://CRAN.R-project.org/package=splines2)
[![Build Status](https://travis-ci.org/wenjie2wang/splines2.svg?branch=master)](https://travis-ci.org/wenjie2wang/splines2)

You can install the released version from [CRAN](https://CRAN.R-project.org/package=splines2).

```R
install.packages("splines2", dependencies = TRUE)
```


## Development

[![Build Status](https://travis-ci.org/wenjie2wang/splines2.svg?branch=dev)](https://travis-ci.org/wenjie2wang/splines2)

The latest version of package is under development
at [GitHub](https://github.com/wenjie2wang/splines2) in branch 'dev'.  You may
consider installing the latest version with the help of **devtools** if it is
able to pass the building check by Travis CI.

```R
if (! require(devtools)) install.packages("devtools", dependencies = TRUE)
devtools::install_git("git://github.com/wenjie2wang/splines2.git", branch = "dev")
```


## Get Started

- [Package vignettes](http://wenjie-stat.me/splines2/)
  provides a quick demonstration for the basic usage of main functions.

- [Package help manual](http://wenjie-stat.me/splines2/splines2.pdf)
  is also available for more technical details.


## License

The R package **splines2** is free software: You can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or any later version
(at your option).  See
the [GNU General Public License](http://www.gnu.org/licenses/) for details.

The R package **splines2** is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.
