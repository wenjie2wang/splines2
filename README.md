splines2
================

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/splines2)](https://CRAN.R-project.org/package=splines2)
[![Downloads from the RStudio CRAN
mirror](https://cranlogs.r-pkg.org/badges/splines2)](https://CRAN.R-project.org/package=splines2)
[![Build
Status](https://travis-ci.org/wenjie2wang/splines2.svg?branch=master)](https://travis-ci.org/wenjie2wang/splines2)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/bvoso7nxchg1incb/branch/master?svg=true)](https://ci.appveyor.com/project/wenjie2wang/splines2)
[![codecov](https://codecov.io/gh/wenjie2wang/splines2/branch/master/graph/badge.svg)](https://codecov.io/gh/wenjie2wang/splines2)

The R package **splines2** is a supplementary package on splines
providing functions constructing B-splines, integral of B-splines,
monotone splines (M-splines) and its integral (I-splines), convex
splines (C-splines), and their derivatives of given order.

## Installation of CRAN Version

You can install the released version from
[CRAN](https://CRAN.R-project.org/package=splines2).

``` r
install.packages("splines2")
```

## Development

The latest version of package is under development at
[GitHub](https://github.com/wenjie2wang/splines2). If it is able to pass
the building check by Travis CI, one may install it by

``` r
if (! require(remotes)) install.packages("remotes")
remotes::install_github("wenjie2wang/splines2")
```

## Getting Started

  - [Online documentation](https://wenjie-stat.me/splines2) provides
    reference for all functions
  - Package vignette for a quick
    [demonstration](https://wenjie-stat.me/splines2/articles/splines2-intro)
    for the common usage through examples

## Performance

Since v0.3.0, the implementation of the main functions has been
rewritten in C++ with the help of the **Rcpp** and **RcppArmadillo**
package. The computational performance has thus been boosted.

Some benchmarks with the **splines** package (version 4.0.1) are
provided for reference as follows:

``` r
library(microbenchmark)
library(splines)
library(splines2)

x <- seq.int(0, 1, 0.001)
degree <- 3
ord <- degree + 1
knots <- seq.int(0.1, 0.9, 0.1)
b_knots <- c(0, 1)
all_knots <- sort(c(knots, rep(range(x), degree + 1)))

## check equivalency of outputs
my_check <- function(values) {
    all(sapply(values[- 1], function(x) {
        all.equal(unclass(values[[1]]), x, check.attributes = FALSE)
    }))
}
```

For B-splines, function `splines2::bSpline()` provides equivalent
results with `splines::bs()` and `splines::splineDesign()`, and is about
3x faster than `bs()` and 2x faster than `splineDesign()`.

``` r
## B-splines
microbenchmark(
    "splines::bs" = bs(x, knots = knots, degree = degree,
                       intercept = TRUE, Boundary.knots = b_knots),
    "splines::splineDesign" = splineDesign(x, knots = all_knots, ord = ord),
    "splines2::bSpline" = bSpline(x, knots = knots, degree = degree,
                                  intercept = TRUE, Boundary.knots = b_knots),
    check = my_check,
    times = 1e3
)
```

``` 
Unit: microseconds
                  expr     min      lq   mean median     uq    max neval cld
           splines::bs 334.774 362.377 418.69 376.99 426.49 3338.6  1000   c
 splines::splineDesign 204.556 215.318 255.04 221.62 242.95 3228.8  1000  b 
     splines2::bSpline  84.762  92.908 118.51  97.21 110.27 2708.6  1000 a  
```

Similarly, for derivatives of B-splines, `splines2::dbs()` provides
equivalent results with `splines::splineDesign()`, and is more than 2x
faster.

``` r
## Derivatives of B-splines
derivs <- 2
microbenchmark(
    "splines::splineDesign" = splineDesign(x, knots = all_knots,
                                           ord = ord, derivs = derivs),
    "splines2::dbs" = dbs(x, derivs = derivs, knots = knots, degree = degree,
                          intercept = TRUE, Boundary.knots = b_knots),
    check = my_check,
    times = 1e3
)
```

    Unit: microseconds
                      expr     min      lq   mean  median     uq    max neval cld
     splines::splineDesign 275.449 285.517 324.67 292.159 309.56 5843.7  1000   b
             splines2::dbs  88.451  94.801 127.73  98.968 106.42 2439.0  1000  a 

The **splines** package does not provide function producing integrals of
B-splines. So we instead performed a comparison with package **ibs**
(version 1.4), where the function `ibs::ibs()` was also implemented in
**Rcpp**.

``` r
## integrals of B-splines
set.seed(123)
coef_sp <- rnorm(length(all_knots) - ord)
microbenchmark(
    "ibs::ibs" = ibs::ibs(x, knots = all_knots, ord = ord, coef = coef_sp),
    "splines2::ibs" =
        as.numeric(splines2::ibs(x, knots = knots, degree = degree,
                                 intercept = TRUE, Boundary.knots = b_knots) %*%
                   coef_sp)
    ,
    check = my_check,
    times = 1e3
)
```

    Unit: microseconds
              expr     min      lq   mean  median      uq      max neval cld
          ibs::ibs 2417.99 2637.78 3216.3 3183.48 3307.88 111745.9  1000   b
     splines2::ibs  264.33  314.82  391.3  337.44  356.06   3577.3  1000  a 

The function `ibs::ibs()` returns the integrated B-splines instead of
the integrals of spline bases. So we applied the same coefficients to
the bases from `splines2::ibs()` for equivalent results, which was still
much faster than `ibs::ibs()`.

## License

[GNU General Public License (v3)](https://www.gnu.org/licenses/)
