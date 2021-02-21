splines2
================

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/splines2)](https://CRAN.R-project.org/package=splines2)
[![Downloads from the RStudio CRAN
mirror](https://cranlogs.r-pkg.org/badges/splines2)](https://CRAN.R-project.org/package=splines2)
[![Build
Status](https://github.com/wenjie2wang/splines2/workflows/R-CMD-check/badge.svg)](https://github.com/wenjie2wang/splines2/actions)
[![codecov](https://codecov.io/gh/wenjie2wang/splines2/branch/main/graph/badge.svg)](https://codecov.io/gh/wenjie2wang/splines2)

The R package **splines2** (version 0.4.2) provides functions to
construct basis matrix of

-   B-splines
-   M-splines
-   I-splines
-   convex splines (C-splines)
-   periodic M-splines
-   natural cubic splines
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
remotes::install_github("wenjie2wang/splines2", upgrade = "never")
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

Some benchmarks are provided for reference as follows:

``` r
library(microbenchmark)
library(splines)
library(splines2)

x <- seq.int(0, 1, 1e-3)
degree <- 3
ord <- degree + 1
internal_knots <- seq.int(0.1, 0.9, 0.1)
boundary_knots <- range(x)
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
3x faster than `bs()` and 2x faster than `splineDesign()`.

``` r
## B-splines
microbenchmark(
    "splines::bs" = bs(x, knots = internal_knots, degree = degree,
                       intercept = TRUE),
    "splines::splineDesign" = splineDesign(x, knots = all_knots, ord = ord),
    "splines2::bSpline" = bSpline(x, knots = internal_knots, degree = degree,
                                  intercept = TRUE),
    check = my_check,
    times = 1e3
)
```

    Unit: microseconds
                      expr     min      lq   mean  median      uq    max neval cld
               splines::bs 333.282 341.558 366.29 346.553 354.520 2357.2  1000   c
     splines::splineDesign 204.503 207.562 234.71 208.990 212.632 2306.0  1000  b 
         splines2::bSpline  84.637  89.913 106.48  93.359  95.321 2171.3  1000 a  

Similarly, for derivatives of B-splines, `splines2::dbs()` provides
equivalent results with `splines::splineDesign()`, and is more than 2x
faster.

``` r
## Derivatives of B-splines
derivs <- 2
microbenchmark(
    "splines::splineDesign" = splineDesign(x, knots = all_knots,
                                           ord = ord, derivs = derivs),
    "splines2::dbs" = dbs(x, derivs = derivs, knots = internal_knots,
                          degree = degree, intercept = TRUE),
    check = my_check,
    times = 1e3
)
```

    Unit: microseconds
                      expr     min      lq   mean  median     uq      max neval cld
     splines::splineDesign 273.854 276.226 302.43 277.329 283.99   2598.4  1000   a
             splines2::dbs  94.245  96.678 261.68  99.937 102.39 139803.9  1000   a

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
        splines2::ibs(x, knots = internal_knots, degree = degree,
                      intercept = TRUE) %*% coef_sp
    ),
    check = my_check,
    times = 1e3
)
```

    Unit: microseconds
              expr     min      lq   mean  median      uq     max neval cld
          ibs::ibs 2492.85 2556.29 3057.2 3180.58 3251.46 17812.4  1000   b
     splines2::ibs  266.24  308.96  336.3  330.91  338.48  3594.9  1000  a 

The function `ibs::ibs()` returns the integrated B-splines instead of
the integrals of spline basis functions. So we applied the same
coefficients to the basis functions from `splines2::ibs()` for
equivalent results, which was still much faster than `ibs::ibs()`.

For natural cubic splines (based on B-splines), `splines::ns()` uses QR
decomposition to find the null space of the second derivatives of
B-spline basis functions at boundary knots, while
`splines2::naturalSpline()` utilizes the close-form null space derived
from the second derivatives of cubic B-splines, which produces
nonnegative basis functions (within boundary) and is more
computationally efficient.

``` r
microbenchmark(
    "splines::ns" = ns(x, knots = internal_knots, intercept = TRUE),
    "splines2::naturalSpline" = naturalSpline(
        x, knots = internal_knots, intercept = TRUE
    ),
    times = 1e3
)
```

    Unit: microseconds
                        expr    min     lq   mean median     uq     max neval cld
                 splines::ns 616.82 627.31 728.50  633.5 648.72 16117.6  1000   b
     splines2::naturalSpline 119.77 128.44 152.71  139.8 143.37  2734.5  1000  a 

The function `mSpline()` produces periodic spline basis functions (based
on M-splines) when `periodic = TRUE` is specified. The
`splines::periodicSpline()` returns a periodic interpolation spline
(based on B-splines) instead of basis matrix. So we performed a
comparison with package **pbs** (version `r packageVersion("pbs")`),
where the function `pbs::pbs()` produces a basis matrix of periodic
B-spline by using `splines::spline.des()` (a wrapper function of
`splines::splineDesign()`).

``` r
microbenchmark(
    "pbs::pbs" = pbs::pbs(x, knots = internal_knots, degree = degree,
                          intercept = TRUE, periodic = TRUE),
    "splines2::mSpline" = mSpline(
        x, knots = internal_knots, degree = degree,
        intercept = TRUE, periodic = TRUE
    ),
    times = 1e3
)
```

    Unit: microseconds
                  expr    min     lq   mean median     uq    max neval cld
              pbs::pbs 422.99 431.45 490.94 435.88 443.24 3660.3  1000   b
     splines2::mSpline 121.80 126.09 147.83 135.34 138.27 2892.8  1000  a 

<details>
<summary>
Session Information for Benchmarks
</summary>

``` r
sessionInfo()
```

    R version 4.0.3 (2020-10-10)
    Platform: x86_64-pc-linux-gnu (64-bit)
    Running under: Arch Linux

    Matrix products: default
    BLAS:   /usr/lib/libopenblasp-r0.3.13.so
    LAPACK: /usr/lib/liblapack.so.3.9.0

    locale:
     [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
     [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
     [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
    [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

    attached base packages:
    [1] splines   stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
    [1] splines2_0.4.2       microbenchmark_1.4-7

    loaded via a namespace (and not attached):
     [1] Rcpp_1.0.6        knitr_1.31        magrittr_2.0.1    MASS_7.3-53.1     ibs_1.4          
     [6] lattice_0.20-41   rlang_0.4.10      multcomp_1.4-16   stringr_1.4.0     tools_4.0.3      
    [11] grid_4.0.3        xfun_0.21         TH.data_1.0-10    htmltools_0.5.1.1 yaml_2.2.1       
    [16] survival_3.2-7    digest_0.6.27     Matrix_1.3-2      codetools_0.2-16  evaluate_0.14    
    [21] rmarkdown_2.6     sandwich_3.0-0    stringi_1.5.3     compiler_4.0.3    pbs_1.1          
    [26] mvtnorm_1.1-1     zoo_1.8-8        

</details>

## License

[GNU General Public License](https://www.gnu.org/licenses/) (≥ 3)
