# splines2

The R package **splines2** is a complementary package on splines providing
functions constructing B-splines, monotone splines (M-splines) and its integral
(I-splines), convex splines (C-splines), and integral of B-splines. Piecewise
constant basis of degree zero is allowed for B-spline and M-spline basis.


## Installation of CRAN Version

[![CRAN_Status_Badge][1]][3]
[![Build Status][4]][5]
[![codecov][codecov-master]][codecov]
[![Downloads from the RStudio CRAN mirror][2]][3]

You can install the released version from [CRAN][3].

```R
install.packages("splines2")
```


## Development

[![Build Status][6]][5]
[![codecov][codecov-dev]](codecov)


The latest version of package is under development at [GitHub][7] in branch
'dev'.  If it is able to pass the building check by Travis CI, you may consider
installing it with the help of **devtools** by

```R
devtools::install_github("wenjie2wang/splines2", ref = "dev")
```

or cloning this reposotory to local and install by makefile as follows:

```
git clone https://github.com/wenjie2wang/splines2.git
cd splines2
make install
```


## Get Started

- [Package vignettes][8]
  provides a quick demonstration for the basic usage of main functions.

- [Package help manual][9] is also available for more technical details.


## License

The R package **splines2** is free software: You can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or any later version
(at your option).  See the [GNU General Public License][10] for details.

The R package **splines2** is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.


[1]: http://www.r-pkg.org/badges/version/splines2
[2]: http://cranlogs.r-pkg.org/badges/splines2
[3]: https://CRAN.R-project.org/package=splines2
[4]: https://travis-ci.org/wenjie2wang/splines2.svg?branch=master
[5]: https://travis-ci.org/wenjie2wang/splines2
[6]: https://travis-ci.org/wenjie2wang/splines2.svg?branch=dev
[7]: https://github.com/wenjie2wang/splines2
[8]: http://wenjie-stat.me/splines2/
[9]: http://wenjie-stat.me/splines2/splines2.pdf
[10]: http://www.gnu.org/licenses/
[codecov]: https://codecov.io/gh/wenjie2wang/splines2
[codecov-master]: https://codecov.io/gh/wenjie2wang/splines2/branch/master/graph/badge.svg
[codecov-dev]: https://codecov.io/gh/wenjie2wang/splines2/branch/dev/graph/badge.svg
