splines2
================

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/splines2)](https://CRAN.R-project.org/package=splines2)
[![Total_Downloads](https://cranlogs.r-pkg.org/badges/grand-total/splines2)](https://CRAN.R-project.org/package=splines2)
[![Downloads from the RStudio CRAN
mirror](https://cranlogs.r-pkg.org/badges/splines2)](https://CRAN.R-project.org/package=splines2)
[![Build
Status](https://github.com/wenjie2wang/splines2/workflows/R-CMD-check/badge.svg)](https://github.com/wenjie2wang/splines2/actions)
[![codecov](https://codecov.io/gh/wenjie2wang/splines2/branch/main/graph/badge.svg)](https://codecov.io/gh/wenjie2wang/splines2)
[![JDS](https://img.shields.io/badge/JDS-10.6339%2F21--JDS1020-brightgreen)](https://doi.org/10.6339/21-JDS1020)

The R package **splines2** is intended to be a user-friendly
*supplement* to the base package **splines**.

## Features

The package **splines2** (version 0.4.5) provides functions to construct
basis matrices of

-   B-splines
-   M-splines
-   I-splines
-   convex splines (C-splines)
-   periodic M-splines
-   natural cubic splines
-   generalized Bernstein polynomials
-   their integrals (except C-splines) and derivatives of given order by
    closed-form recursive formulas

In addition to the R interface, **splines2** provides a C++ header-only
library integrated with **Rcpp**, which allows the construction of
spline basis functions directly in C++ with the help of **Rcpp** and
**RcppArmadillo**. Thus, it can also be treated as one of the **Rcpp\***
packages. A toy example package that uses the C++ interface is available
[here](https://github.com/wenjie2wang/example-pkg-Rcpp-splines2).

## Installation of CRAN Version

You can install the released version from
[CRAN](https://CRAN.R-project.org/package=splines2).

``` r
install.packages("splines2")
```

## Development

The latest version of the package is under development at
[GitHub](https://github.com/wenjie2wang/splines2). If it is able to pass
the automated package checks, one may install it by

``` r
if (! require(remotes)) install.packages("remotes")
remotes::install_github("wenjie2wang/splines2", upgrade = "never")
```

## Getting Started

The [Online document](https://wwenjie.org/splines2) provides a reference
for all functions and contains the following vignettes:

-   [Demonstration of the common usages in R through
    examples](https://wwenjie.org/splines2/articles/splines2-intro).
-   [Introduction to the usage with
    Rcpp](https://wwenjie.org/splines2/articles/splines2-wi-rcpp)

## Performance

Since v0.3.0, the implementation of the main functions has been
rewritten in C++ with the help of the **Rcpp** and **RcppArmadillo**
packages. The computational performance has thus been boosted and
comparable with the function `splines::splineDesign()`.

Some quick micro-benchmarks are provided for reference as follows:

``` r
library(microbenchmark)
library(splines)
library(splines2)

set.seed(123)
x <- runif(1e3)
degree <- 3
ord <- degree + 1
internal_knots <- seq.int(0.1, 0.9, 0.1)
boundary_knots <- c(0, 1)
all_knots <- sort(c(internal_knots, rep(boundary_knots, ord)))

## check equivalency of outputs
my_check <- function(values) {
    all(sapply(values[- 1], function(x) {
        all.equal(unclass(values[[1]]), x, check.attributes = FALSE)
    }))
}
```

For B-splines, function `splines2::bSpline()` provides equivalent
results with `splines::bs()` and `splines::splineDesign()`, and is about
3x faster than `bs()` and 2x faster than `splineDesign()` for this
example.

``` r
## B-splines
microbenchmark(
    "splines::bs" = bs(x, knots = internal_knots, degree = degree,
                       intercept = TRUE, Boundary.knots = boundary_knots),
    "splines::splineDesign" = splineDesign(x, knots = all_knots, ord = ord),
    "splines2::bSpline" = bSpline(
        x, knots = internal_knots, degree = degree,
        intercept = TRUE, Boundary.knots = boundary_knots
    ),
    check = my_check,
    times = 1e3,
    unit = "relative"
)
```

    Unit: relative
                      expr    min     lq   mean median     uq    max neval cld
               splines::bs 3.6436 3.5392 3.5667 3.4591 3.4965 1.2118  1000   c
     splines::splineDesign 2.2156 2.1265 2.1414 2.0520 2.1235 1.0717  1000  b 
         splines2::bSpline 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000  1000 a  

Similarly, for derivatives of B-splines, `splines2::dbs()` provides
equivalent results with `splines::splineDesign()`, and is about 2x
faster.

``` r
## Derivatives of B-splines
derivs <- 2
microbenchmark(
    "splines::splineDesign" = splineDesign(x, knots = all_knots,
                                           ord = ord, derivs = derivs),
    "splines2::dbs" = dbs(x, derivs = derivs, knots = internal_knots,
                          degree = degree, intercept = TRUE,
                          Boundary.knots = boundary_knots),
    check = my_check,
    times = 1e3,
    unit = "relative"
)
```

    Unit: relative
                      expr    min     lq   mean median     uq    max neval cld
     splines::splineDesign 2.6213 2.4951 2.4608 2.4362 2.4634 1.1087  1000   b
             splines2::dbs 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000  1000  a 

The **splines** package does not contain an implementation for integrals
of B-splines. Thus, we performed a comparison with package **ibs**
(version `r packageVersion("ibs")`), where the function `ibs::ibs()` was
also implemented in **Rcpp**.

``` r
## integrals of B-splines
set.seed(123)
coef_sp <- rnorm(length(all_knots) - ord)
microbenchmark(
    "ibs::ibs" = ibs::ibs(x, knots = all_knots, ord = ord, coef = coef_sp),
    "splines2::ibs" = as.numeric(
        splines2::ibs(x, knots = internal_knots, degree = degree,
                      intercept = TRUE, Boundary.knots = boundary_knots) %*%
        coef_sp
    ),
    check = my_check,
    times = 1e3,
    unit = "relative"
)
```

    Unit: relative
              expr    min     lq   mean median     uq    max neval cld
          ibs::ibs 9.2171 8.2285 6.5033 9.5655 9.4491 18.632  1000   b
     splines2::ibs 1.0000 1.0000 1.0000 1.0000 1.0000  1.000  1000  a 

The function `ibs::ibs()` returns the integrated B-splines instead of
the integrals of spline basis functions. Thus, we applied the same
coefficients to the basis functions from `splines2::ibs()` for
equivalent results, which was still much faster than `ibs::ibs()`.

For natural cubic splines (based on B-splines), `splines::ns()` uses the
QR decomposition to find the null space of the second derivatives of
B-spline basis functions at boundary knots, while
`splines2::naturalSpline()` utilizes the closed-form null space derived
from the second derivatives of cubic B-splines, which produces
nonnegative basis functions (within boundary) and is more
computationally efficient.

``` r
microbenchmark(
    "splines::ns" = ns(x, knots = internal_knots, intercept = TRUE,
                       Boundary.knots = boundary_knots),
    "splines2::naturalSpline" = naturalSpline(
        x, knots = internal_knots, intercept = TRUE,
        Boundary.knots = boundary_knots
    ),
    times = 1e3,
    unit = "relative"
)
```

    Unit: relative
                        expr    min     lq   mean median     uq    max neval cld
                 splines::ns 5.1898 4.9857 4.7973  4.736 4.5995 1.6821  1000   b
     splines2::naturalSpline 1.0000 1.0000 1.0000  1.000 1.0000 1.0000  1000  a 

The function `mSpline()` produces periodic spline basis functions (based
on M-splines) when `periodic = TRUE` is specified. The
`splines::periodicSpline()` returns a periodic interpolation spline
(based on B-splines) instead of basis matrix. Thus, we performed a
comparison with package **pbs** (version `r packageVersion("pbs")`),
where the function `pbs::pbs()` produces a basis matrix of periodic
B-spline by using `splines::spline.des()` (a wrapper function of
`splines::splineDesign()`).

``` r
microbenchmark(
    "pbs::pbs" = pbs::pbs(x, knots = internal_knots, degree = degree,
                          intercept = TRUE, periodic = TRUE,
                          Boundary.knots = boundary_knots),
    "splines2::mSpline" = mSpline(
        x, knots = internal_knots, degree = degree, intercept = TRUE,
        Boundary.knots = boundary_knots, periodic = TRUE
    ),
    times = 1e3,
    unit = "relative"
)
```

    Unit: relative
                  expr    min    lq   mean median     uq    max neval cld
              pbs::pbs 3.3766 3.229 3.1478 3.1355 3.1166 1.4485  1000   b
     splines2::mSpline 1.0000 1.000 1.0000 1.0000 1.0000 1.0000  1000  a 

<details>
<summary>
Session Information for Benchmarks
</summary>

``` r
sessionInfo()
```

    R version 4.1.0 (2021-05-18)
    Platform: x86_64-pc-linux-gnu (64-bit)
    Running under: Arch Linux

    Matrix products: default
    BLAS:   /usr/lib/libopenblasp-r0.3.17.so
    LAPACK: /usr/lib/liblapack.so.3.10.0

    locale:
     [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
     [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
     [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
    [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

    attached base packages:
    [1] splines   stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
    [1] splines2_0.4.5       microbenchmark_1.4-7

    loaded via a namespace (and not attached):
     [1] Rcpp_1.0.7        mvtnorm_1.1-2     lattice_0.20-44   codetools_0.2-18  ibs_1.4          
     [6] zoo_1.8-9         digest_0.6.27     MASS_7.3-54       grid_4.1.0        magrittr_2.0.1   
    [11] evaluate_0.14     rlang_0.4.11      stringi_1.7.3     multcomp_1.4-17   Matrix_1.3-4     
    [16] sandwich_3.0-1    rmarkdown_2.10    TH.data_1.0-10    tools_4.1.0       stringr_1.4.0    
    [21] survival_3.2-11   xfun_0.25         yaml_2.2.1        compiler_4.1.0    pbs_1.1          
    [26] htmltools_0.5.1.1 knitr_1.33       

</details>

## License

[GNU General Public License](https://www.gnu.org/licenses/) (≥ 3)
