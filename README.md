splines2
================

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/splines2)](https://CRAN.R-project.org/package=splines2)
[![Downloads from the RStudio CRAN
mirror](https://cranlogs.r-pkg.org/badges/splines2)](https://CRAN.R-project.org/package=splines2)
[![Build
Status](https://github.com/wenjie2wang/splines2/workflows/R-CMD-check/badge.svg)](https://github.com/wenjie2wang/splines2/actions)
[![codecov](https://codecov.io/gh/wenjie2wang/splines2/branch/main/graph/badge.svg)](https://codecov.io/gh/wenjie2wang/splines2)

The R package **splines2** provides functions to construct basis matrix
of

-   B-splines
-   M-splines
-   I-splines
-   convex splines (C-splines)
-   generalized Bernstein polynomials
-   their integrals (except C-splines) and derivatives of given order by
    close-form recursive formulas

In addition to the R interface, **splines2** also provides a C++
header-only library integrated with **Rcpp**, which allows construction
of spline basis matrix directly in C++ with the help of **Rcpp** and
**RcppArmadillo**. So it can also be treated as one of the **Rcpp\***
packages. A toy example package that uses the C++ interface is available
[here](https://github.com/wenjie2wang/example-pkg-Rcpp-splines2).

## Installation of CRAN Version

You can install the released version from
[CRAN](https://CRAN.R-project.org/package=splines2).

``` r
install.packages("splines2")
```

## Development

The latest version of package is under development at
[GitHub](https://github.com/wenjie2wang/splines2). If it is able to pass
the automated package checks, one may install it by

``` r
if (! require(remotes)) install.packages("remotes")
remotes::install_github("wenjie2wang/splines2")
```

## Getting Started

[Online document](https://wwenjie.org/splines2) provides reference for
all functions and contains the following vignettes:

-   [Demonstration of the common usages in R through
    examples](https://wwenjie.org/splines2/articles/splines2-intro).
-   [Introduction to the usage with
    Rcpp](https://wwenjie.org/splines2/articles/splines2-wi-rcpp)

## Performance

Since v0.3.0, the implementation of the main functions has been
rewritten in C++ with the help of the **Rcpp** and **RcppArmadillo**
package. The computational performance has thus been boosted.

Some benchmarks with the **splines** package (version 4.0.3) are
provided for reference as follows:

``` r
library(microbenchmark)
library(splines)
library(splines2)

x <- seq.int(0, 1, 0.001)
degree <- 3
ord <- degree + 1
knots <- seq.int(0.1, 0.9, 0.1)
b_knots <- range(x)
all_knots <- sort(c(knots, rep(b_knots, ord)))

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

    Unit: microseconds
                      expr     min      lq   mean median      uq    max neval cld
               splines::bs 333.624 345.985 381.45 353.56 369.332 2459.1  1000   c
     splines::splineDesign 204.708 208.529 239.15 211.58 220.958 2328.3  1000  b 
         splines2::bSpline  84.764  90.584 101.50  93.84  97.044 2113.6  1000 a  

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
                      expr     min      lq   mean  median      uq    max neval cld
     splines::splineDesign 273.311 277.778 324.64 283.183 292.754 5308.8  1000   b
             splines2::dbs  87.589  92.416 103.81  96.128  99.986 2355.2  1000  a 

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
    "splines2::ibs" = as.numeric(
        splines2::ibs(x, knots = knots, degree = degree,
                      intercept = TRUE, Boundary.knots = b_knots) %*% coef_sp
    ),
    check = my_check,
    times = 1e3
)
```

    Unit: microseconds
              expr     min      lq    mean  median     uq      max neval cld
          ibs::ibs 2361.66 2601.33 3188.16 3115.31 3337.0 114069.3  1000   b
     splines2::ibs  251.87  309.05  377.52  333.76  345.6   3806.2  1000  a 

The function `ibs::ibs()` returns the integrated B-splines instead of
the integrals of spline bases. So we applied the same coefficients to
the bases from `splines2::ibs()` for equivalent results, which was still
much faster than `ibs::ibs()`.

## License

[GNU General Public License](https://www.gnu.org/licenses/) (≥ 3)
