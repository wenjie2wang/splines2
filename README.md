splines2
================

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/splines2)](https://CRAN.R-project.org/package=splines2)
[![Total_Downloads](https://cranlogs.r-pkg.org/badges/grand-total/splines2)](https://CRAN.R-project.org/package=splines2)
[![Downloads from the RStudio CRAN
mirror](https://cranlogs.r-pkg.org/badges/splines2)](https://CRAN.R-project.org/package=splines2)
[![Build
Status](https://github.com/wenjie2wang/splines2/workflows/R-CMD-check/badge.svg)](https://github.com/wenjie2wang/splines2/actions)
[![codecov](https://codecov.io/gh/wenjie2wang/splines2/branch/main/graph/badge.svg)](https://app.codecov.io/gh/wenjie2wang/splines2)
[![JDS](https://img.shields.io/badge/JDS-10.6339%2F21--JDS1020-brightgreen)](https://doi.org/10.6339/21-JDS1020)

**Package website**: [release](https://wwenjie.org/splines2) \|
[development](https://wwenjie.org/splines2/dev)

The R package **splines2** is intended to be a user-friendly
*supplementary* package to the base package **splines**.

## Features

The package **splines2** provides functions to construct basis matrices
of

- B-splines
- M-splines
- I-splines
- convex splines (C-splines)
- periodic splines
- natural cubic splines
- generalized Bernstein polynomials
- their integrals (except C-splines) and derivatives of given order by
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

- [Demonstration of the common usages in R through
  examples](https://wwenjie.org/splines2/articles/splines2-intro).
- [Introduction to the usage with
  Rcpp](https://wwenjie.org/splines2/articles/splines2-wi-rcpp)

## Performance

Since version 0.3.0, the implementation of the main functions has been
rewritten in C++ with the help of the **Rcpp** and **RcppArmadillo**
packages. The computational performance has thus been boosted and
comparable with the function `splines::splineDesign()`.

Some quick micro-benchmarks are provided below for reference.

``` r
library(microbenchmark)
options(microbenchmark.unit="relative")
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
        all.equal(unclass(values[[1]]), unclass(x), check.attributes = FALSE)
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
    check = my_check
)
```

    Unit: relative
                      expr    min     lq   mean median     uq    max neval
               splines::bs 3.5101 2.7677 2.4657 2.8721 2.2682 1.5334   100
     splines::splineDesign 2.4324 1.9049 1.9833 2.1624 1.7941 8.4094   100
         splines2::bSpline 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000   100

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
    check = my_check
)
```

    Unit: relative
                      expr    min     lq   mean median     uq    max neval
     splines::splineDesign 2.7082 2.4949 2.0586 2.0985 2.1164 1.5585   100
             splines2::dbs 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000   100

The **splines** package does not contain an implementation for integrals
of B-splines. Thus, we performed a comparison with package **ibs**
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
    check = my_check
)
```

    Unit: relative
              expr    min     lq   mean median     uq    max neval
          ibs::ibs 53.484 45.639 42.991 42.464 39.945 30.853   100
     splines2::ibs  1.000  1.000  1.000  1.000  1.000  1.000   100

The function `ibs::ibs()` returns the integrated B-splines instead of
the integrals of spline basis functions. Thus, we applied the same
coefficients to the basis functions from `splines2::ibs()` for
equivalent results, which was still much faster than `ibs::ibs()`.

For natural cubic splines (based on B-splines), `splines::ns()` uses the
QR decomposition to find the null space of the second derivatives of
B-spline basis functions at boundary knots, while `splines2::nsp()`
utilizes the closed-form null space derived from the second derivatives
of cubic B-splines, which produces nonnegative basis functions (within
boundary) and is more computationally efficient.

``` r
microbenchmark(
    "splines::ns" = ns(x, knots = internal_knots, intercept = TRUE,
                       Boundary.knots = boundary_knots),
    "splines2::nsp" = nsp(
        x, knots = internal_knots, intercept = TRUE,
        Boundary.knots = boundary_knots
    )
)
```

    Unit: relative
              expr    min     lq   mean median     uq    max neval
       splines::ns 5.8949 5.4165 4.8466 4.5983 5.1892 3.7201   100
     splines2::nsp 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000   100

The functions `bSpline()` and `mSpline()` produce periodic spline basis
functions based on B-splines and M-splines, respectively, when
`periodic = TRUE` is specified. The `splines::periodicSpline()` returns
a periodic interpolation spline (based on B-splines) instead of basis
matrix. We performed a comparison with package **pbs** (version 1.1),
where the function `pbs::pbs()` produces a basis matrix of periodic
B-spline by using `splines::spline.des()`.

``` r
microbenchmark(
    "pbs::pbs" = pbs::pbs(x, knots = internal_knots, degree = degree,
                          intercept = TRUE, periodic = TRUE,
                          Boundary.knots = boundary_knots),
    "splines2::bSpline" = bSpline(
        x, knots = internal_knots, degree = degree, intercept = TRUE,
        Boundary.knots = boundary_knots, periodic = TRUE
    ),
    "splines2::mSpline" = mSpline(
        x, knots = internal_knots, degree = degree, intercept = TRUE,
        Boundary.knots = boundary_knots, periodic = TRUE
    )
)
```

    Unit: relative
                  expr    min     lq   mean median     uq     max neval
              pbs::pbs 4.2867 3.8656 2.6849 3.6822 2.8756 1.00982   100
     splines2::bSpline 1.0000 1.0000 1.0000 1.0000 1.0000 1.00000   100
     splines2::mSpline 1.0787 1.0477 0.8013 1.1933 1.0466 0.11527   100

<details>
<summary>
Session Information for Benchmarks
</summary>

``` r
sessionInfo()
```

    R version 4.4.2 (2024-10-31)
    Platform: x86_64-pc-linux-gnu
    Running under: Arch Linux

    Matrix products: default
    BLAS/LAPACK: /usr/lib/libopenblas.so.0.3;  LAPACK version 3.12.0

    locale:
     [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
     [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
     [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
    [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

    time zone: America/Indianapolis
    tzcode source: system (glibc)

    attached base packages:
    [1] splines   stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
    [1] splines2_0.5.4       microbenchmark_1.5.0

    loaded via a namespace (and not attached):
     [1] digest_0.6.37     codetools_0.2-20  ibs_1.4           fastmap_1.2.0     xfun_0.51        
     [6] pbs_1.1           knitr_1.49        htmltools_0.5.8.1 rmarkdown_2.29    cli_3.6.4        
    [11] compiler_4.4.2    tools_4.4.2       evaluate_1.0.3    Rcpp_1.0.14       yaml_2.3.10      
    [16] rlang_1.1.5      

</details>

## License

[GNU General Public License](https://www.gnu.org/licenses/) (≥ 3)
