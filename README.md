splines2
================

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/splines2)](https://CRAN.R-project.org/package=splines2)
[![Downloads from the RStudio CRAN
mirror](https://cranlogs.r-pkg.org/badges/splines2)](https://CRAN.R-project.org/package=splines2)
[![Build
Status](https://github.com/wenjie2wang/splines2/workflows/R-CMD-check/badge.svg)](https://github.com/wenjie2wang/splines2/actions)
[![codecov](https://codecov.io/gh/wenjie2wang/splines2/branch/main/graph/badge.svg)](https://codecov.io/gh/wenjie2wang/splines2)

The R package **splines2** (version 0.4.0.9000) provides functions to
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
library(splines); packageVersion("splines")
```

    [1] '4.0.3'

``` r
library(splines2); packageVersion("splines2")
```

    [1] '0.4.0.9000'

``` r
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
                      expr     min      lq   mean  median      uq    max neval cld
               splines::bs 335.683 344.499 382.18 352.441 368.265 2531.3  1000   c
     splines::splineDesign 206.776 210.522 233.81 213.150 221.888 2340.3  1000  b 
         splines2::bSpline  85.427  91.076 111.01  94.108  96.919 2260.3  1000 a  

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
                      expr     min      lq   mean median     uq    max neval cld
     splines::splineDesign 276.253 279.967 319.62 285.29 299.21 5047.1  1000   b
             splines2::dbs  94.798  98.338 121.30 101.53 106.29 2571.0  1000  a 

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
              expr     min      lq    mean  median      uq      max neval cld
          ibs::ibs 2366.11 2632.01 3195.60 3116.11 3338.40 115479.6  1000   b
     splines2::ibs  269.92  312.98  339.78  345.56  354.21   1246.3  1000  a 

The function `ibs::ibs()` returns the integrated B-splines instead of
the integrals of spline bases. So we applied the same coefficients to
the bases from `splines2::ibs()` for equivalent results, which was still
much faster than `ibs::ibs()`.

For natural cubic splines (based on B-splines), `splines::ns()` uses QR
decomposition to find the null space of the second derivatives of
B-spline basis at boundary knots, while `splines2::naturalSpline()`
utilizes the close-form null space derived from the second derivatives
of cubic B-splines, which produces nonnegative bases (within boundary)
and is more computationally efficient.

``` r
microbenchmark(
    "splines::ns" = ns(x, knots = knots, intercept = TRUE,
                       Boundary.knots = b_knots),
    "splines2::naturalSpline" = naturalSpline(
        x, knots = knots, intercept = TRUE, Boundary.knots = b_knots
    ),
    times = 1e3
)
```

    Unit: microseconds
                        expr    min     lq   mean median     uq    max neval cld
                 splines::ns 622.45 635.51 716.01 646.20 662.50 3544.9  1000   b
     splines2::naturalSpline 130.89 139.24 166.53 149.63 153.49 2826.6  1000  a 

The function `mSpline()` produces periodic spline basis (based on
M-splines) when `periodic = TRUE` is specified. The
`splines::periodicSpline()` returns a periodic interpolation spline
(based on B-splines) instead of basis matrix. So we performed a
comparison with package **pbs** (version 1.1), where the function
`pbs::pbs()` produces a basis matrix of periodic B-spline by using
`splines::spline.des()` (a wrapper function of
`splines::splineDesign()`).

``` r
microbenchmark(
    "pbs::pbs" = pbs::pbs(x, knots = knots, degree = degree, intercept = TRUE,
                          Boundary.knots = b_knots, periodic = TRUE),
    "splines2::mSpline" = mSpline(
        x, knots = knots, degree = degree, intercept = TRUE,
        Boundary.knots = b_knots, periodic = TRUE
    ),
    times = 1e3
)
```

    Unit: microseconds
                  expr    min     lq   mean median     uq    max neval cld
              pbs::pbs 429.31 439.80 498.89 448.09 460.72 3687.4  1000   b
     splines2::mSpline 122.63 129.48 153.46 138.70 143.17 2963.1  1000  a 

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
    BLAS:   /usr/lib/libopenblasp-r0.3.12.so
    LAPACK: /usr/lib/liblapack.so.3.9.0

    locale:
     [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
     [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
     [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
    [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

    attached base packages:
    [1] splines   stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
    [1] splines2_0.4.0.9000  microbenchmark_1.4-7

    loaded via a namespace (and not attached):
     [1] Rcpp_1.0.5       knitr_1.30       magrittr_2.0.1   MASS_7.3-53      ibs_1.4         
     [6] lattice_0.20-41  rlang_0.4.9      multcomp_1.4-15  stringr_1.4.0    tools_4.0.3     
    [11] grid_4.0.3       xfun_0.19        TH.data_1.0-10   htmltools_0.5.0  yaml_2.2.1      
    [16] survival_3.2-7   digest_0.6.27    Matrix_1.2-18    codetools_0.2-16 evaluate_0.14   
    [21] rmarkdown_2.6    sandwich_3.0-0   stringi_1.5.3    compiler_4.0.3   pbs_1.1         
    [26] mvtnorm_1.1-1    zoo_1.8-8       

</details>

## License

[GNU General Public License](https://www.gnu.org/licenses/) (≥ 3)
