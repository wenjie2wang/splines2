# splines2 0.4.3

## New features

* Enabled extended knot sequence that allows multiplicity of internal knots for
  B-splines, M-splines, I-splines, and C-splines in the C++ interface.
* Added type conversion to `BernsteinPoly` and `PeriodicMSpline` objects to the
  C++ interface.

## Minor changes

* Added testing examples for constructing spline basis functions via the C++
  interface.


# splines2 0.4.2

## New features

* Added `knots()` methods to extract internal knots and boundary knots from a
  given *splines2* object.

## Major changes

* Updated the generation of the knot sequence for period M-splines following
  Piegl and Tiller (1997), which relaxed the previous requirement that
  `length(knots) >= degree` to `length(knots) >= degree - 1`.


# splines2 0.4.1

## New features

* Added function `naturalSpline()` providing implementation of nonnegative
  natural cubic splines.
* Added argument `periodic` to function `mSpline()` for periodic M-splines.
* Added argument `integral` to function `mSpline()` for integrals of M-splines
  or periodic M-splines.
* Added `deriv`, `predict`, and `print` method for `naturalSpline` class object.

## Minor changes

* Updated the `deriv()` method for `mSpline` class object for periodic
  M-splines.


# splines2 0.3.1

## Minor changes

* Modified testing examples for CRAN tests on r-patched-solaris-x86 and
  r-oldrel-macos-x86_64.


# splines2 0.3.0

## New features

* Added function `BersteinPoly()` providing implementation of generalized
  Bernstein polynomials.
* Added C++ interface that can be easily integrated with **Rcpp**.

## Major changes

* Changed most implementations from R to C++ with help of **Rcpp** and
  **RcppArmadillo** to boost the performance.

## Minor changes

* Made piece-wise constant basis functions continuous at right boundary knot for
  consistency with spline basis matrix of non-zero degrees.
* Changed the default value of argument `intercept` in function `iSpline()` and
  `cSpline()` to `TRUE` for a complete set of spline basis functions in
  shape-restricted regression.
* Removed the corresponding M-spline basis from attributes of outputs from
  `iSpline()` and `cSpline()`.
* Removed the corresponding B-spline basis from attributes of outputs from
  `bSpline()`.

## Bug fixes

* Fixed `deriv.mSpline()` method for third derivatives of scaled C-splines.


# splines2 0.2.8

## Bug fixes

* Fixed inconsistency of argument `df` for piecewise constant basis functions
  when `knots = NULL`.

## Minor changes

* Rewrote testing suite for using the **tinytest** package instead of
  **testthat**.


# splines2 0.2.7

## Minor changes

* Updated tests for R development version.


# splines2 0.2.6

## Minor changes

* Added checks for any internal knot incorrectly placed outside of the boundary
  knots and added warnings for users' reference.


# splines2 0.2.5

## Minor changes

* Added more tests and increased code coverage.

## Bug fixes

* Fixed evaluation of derivatives of M-splines for a single value. Thanks Ina
  Jazic for reporting the bug and providing possible fix.
* Fixed `deriv.cSpline()` method for derivatives of order greater than two when
  `scale = TRUE`.


# splines2 0.2.4

## New features

* Added function `dbs()` generating derivative of given order of B-splines. It
  is a similar function with `splines::splineDesign()`. However, it provides a
  more user-friendly interface and more consistent handling on `NA`'s.
* Added `deriv()` methods for derivatives of given order of any existing
  **splines2** object that can be generated currently.

## Major changes

* Added argument `derivs` to function `mSpline()` and `iSpline()` for
  derivatives.
* Changed all the classes of object generated for a better dispatching on
  methods.

## Minor changes

* Added tests for all major functions with the help of package **testthat**.

## Bug fixes

* Fixed the generation of splines without any internal knot.


# splines2 0.2.3

## Bug fixes

* Fixed one-piece constant basis for M-splines.


# splines2 0.2.2

## Bug fixes

* Fixed the NA's handling in all the functions constructing spline basis matrix.


# splines2 0.2.1

## New features

* Added function `bSpline()` generating B-spline basis allowing zero degree or
  piecewise constant basis based on function `bs()` in the **splines** package.
* Introduced function `bSpline()` to allow M-splines of degree zero.
* Added function `cSpline()` constructing convex spline (C-spline) basis.
* Added `predict()` methods for `bSpline2` object and `cSpline` object generated
  by `bSpline()` and `cSpline()`, respectively.
* Added `print()` methods for all **splines2** objects developed so far.

## Major changes

* Improved the function `iSpline()` to construct I-spline basis directly from
  B-spline basis instead of M-spline basis.

## Minor changes

* Updated all CRAN URL to a canonical form suggested.


# splines2 0.1.0

## New features

* The first version of **splines2** providing functions constructing M-spline,
  I-spline, and integral of B-spline basis.
