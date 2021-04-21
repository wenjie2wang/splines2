splines2
================

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/splines2)](https://CRAN.R-project.org/package=splines2)
[![Downloads from the RStudio CRAN
mirror](https://cranlogs.r-pkg.org/badges/splines2)](https://CRAN.R-project.org/package=splines2)
[![Build
Status](https://github.com/wenjie2wang/splines2/workflows/R-CMD-check/badge.svg)](https://github.com/wenjie2wang/splines2/actions)
[![codecov](https://codecov.io/gh/wenjie2wang/splines2/branch/main/graph/badge.svg)](https://codecov.io/gh/wenjie2wang/splines2)

The R package **splines2** (version 0.4.3) provides functions to
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

In addition to the R interface, **splines2** provides a C++ header-only
library integrated with **Rcpp**, which allows construction of spline
basis matrics directly in C++ with the help of **Rcpp** and
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
package. The computational performance has thus been boosted and
comparable with the function `splines::splineDesign()`.

Some quick microbenchmarks are provided for reference as follows:

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
    times = 1e3
)
```

    Unit: microseconds
                      expr     min      lq   mean median     uq    max neval cld
               splines::bs 336.731 349.625 375.79 357.79 372.94 2459.4  1000   c
     splines::splineDesign 207.444 211.532 251.32 213.66 223.17 2452.3  1000  b 
         splines2::bSpline  92.542  99.558 110.78 104.36 107.63 2152.3  1000 a  

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
    times = 1e3
)
```

    Unit: microseconds
                      expr    min     lq   mean median     uq    max neval cld
     splines::splineDesign 276.58 281.51 308.49 284.38 300.67 2783.0  1000   b
             splines2::dbs 108.04 115.40 149.81 120.57 124.98 2380.4  1000  a 

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
                      intercept = TRUE, Boundary.knots = boundary_knots) %*%
        coef_sp
    ),
    check = my_check,
    times = 1e3
)
```

    Unit: microseconds
              expr     min      lq    mean  median      uq      max neval cld
          ibs::ibs 2403.63 2766.26 3363.20 3277.27 3456.95 158210.4  1000   b
     splines2::ibs  293.95  340.09  370.32  377.29  387.72   1231.8  1000  a 

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
    "splines::ns" = ns(x, knots = internal_knots, intercept = TRUE,
                       Boundary.knots = boundary_knots),
    "splines2::naturalSpline" = naturalSpline(
        x, knots = internal_knots, intercept = TRUE,
        Boundary.knots = boundary_knots
    ),
    times = 1e3
)
```

    Unit: microseconds
                        expr    min     lq   mean median     uq    max neval cld
                 splines::ns 628.24 649.58 742.05 663.07 681.58 3486.3  1000   b
     splines2::naturalSpline 126.33 133.88 154.71 143.34 147.69 2677.3  1000  a 

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
                          intercept = TRUE, periodic = TRUE,
                          Boundary.knots = boundary_knots),
    "splines2::mSpline" = mSpline(
        x, knots = internal_knots, degree = degree, intercept = TRUE,
        Boundary.knots = boundary_knots, periodic = TRUE
    ),
    times = 1e3
)
```

    Unit: microseconds
                  expr    min     lq   mean median     uq     max neval cld
              pbs::pbs 428.75 440.79 523.18 449.95 468.11 10239.2  1000   b
     splines2::mSpline 123.40 133.94 150.27 142.59 147.99  2754.9  1000  a 

<details>
<summary>
Session Information for Benchmarks
</summary>

``` r
sessionInfo()
```

    R version 4.0.5 (2021-03-31)
    Platform: x86_64-pc-linux-gnu (64-bit)
    Running under: Arch Linux

    Matrix products: default
    BLAS:   /usr/lib/libopenblasp-r0.3.13.so
    LAPACK: /usr/lib/liblapack.so.3.9.1

    locale:
     [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
     [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
     [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
    [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

    attached base packages:
    [1] splines   stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
    [1] splines2_0.4.3       microbenchmark_1.4-7

    loaded via a namespace (and not attached):
     [1] Rcpp_1.0.6        mvtnorm_1.1-1     lattice_0.20-41   codetools_0.2-18  ibs_1.4          
     [6] zoo_1.8-9         digest_0.6.27     MASS_7.3-53.1     grid_4.0.5        magrittr_2.0.1   
    [11] evaluate_0.14     rlang_0.4.10      stringi_1.5.3     multcomp_1.4-16   Matrix_1.3-2     
    [16] sandwich_3.0-0    rmarkdown_2.7     TH.data_1.0-10    tools_4.0.5       stringr_1.4.0    
    [21] survival_3.2-10   xfun_0.22         yaml_2.2.1        compiler_4.0.5    pbs_1.1          
    [26] htmltools_0.5.1.1 knitr_1.32       

</details>

## License

[GNU General Public License](https://www.gnu.org/licenses/) (≥ 3)
