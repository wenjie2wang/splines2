splines2
================

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/splines2)](https://CRAN.R-project.org/package=splines2)
[![Downloads from the RStudio CRAN
mirror](https://cranlogs.r-pkg.org/badges/splines2)](https://CRAN.R-project.org/package=splines2)
[![Build
Status](https://github.com/wenjie2wang/splines2/workflows/R-CMD-check/badge.svg)](https://github.com/wenjie2wang/splines2/actions)
[![codecov](https://codecov.io/gh/wenjie2wang/splines2/branch/main/graph/badge.svg)](https://codecov.io/gh/wenjie2wang/splines2)

The R package **splines2** (version 0.4.1) provides functions to
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
               splines::bs 336.952 345.183 375.79 350.601 360.875 2505.0  1000   c
     splines::splineDesign 205.412 208.695 234.77 210.150 214.097 2192.4  1000  b 
         splines2::bSpline  85.711  90.792 103.96  94.229  96.755 2143.7  1000 a  

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
     splines::splineDesign 274.138 277.414 315.71 278.56 285.50 2538.4  1000   b
             splines2::dbs  95.243  97.353 115.49 100.96 103.44 2532.8  1000  a 

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
          ibs::ibs 2359.14 2533.36 3173.31 3155.07 3214.27 155132.5  1000   b
     splines2::ibs  265.85  298.67  330.44  333.34  346.81   1267.2  1000  a 

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
                 splines::ns 622.39 632.32 721.41 638.08 649.01 3409.0  1000   b
     splines2::naturalSpline 120.08 129.01 144.23 140.33 143.45 2661.1  1000  a 

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
              pbs::pbs 428.00 437.30 497.21 442.29 452.67 9996.2  1000   b
     splines2::mSpline 122.13 126.25 155.73 135.47 138.67 2934.3  1000  a 

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
    [1] splines2_0.4.1       microbenchmark_1.4-7

    loaded via a namespace (and not attached):
     [1] Rcpp_1.0.5       knitr_1.30       magrittr_2.0.1   MASS_7.3-53      ibs_1.4         
     [6] lattice_0.20-41  rlang_0.4.10     multcomp_1.4-15  stringr_1.4.0    tools_4.0.3     
    [11] grid_4.0.3       xfun_0.19        TH.data_1.0-10   htmltools_0.5.0  yaml_2.2.1      
    [16] survival_3.2-7   digest_0.6.27    Matrix_1.3-0     codetools_0.2-16 evaluate_0.14   
    [21] rmarkdown_2.6    sandwich_3.0-0   stringi_1.5.3    compiler_4.0.3   pbs_1.1         
    [26] mvtnorm_1.1-1    zoo_1.8-8       

</details>

## License

[GNU General Public License](https://www.gnu.org/licenses/) (≥ 3)
