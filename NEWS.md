# splines2 0.5.0.9000

## New features

* Added a new function named `nsk()` for natural cubic spline basis functions
  following the `survival::nsk()` function suggested by Dr. Terry Therneau.
* Added `plot()` methods to quickly visualize the spline basis functions.
* Added `$` method to extract an attribute of the returned `splines2` object.
* Added a new argument named `periodic` to `bSpline()` for periodic B-splines
  and a new class named `PeriodicBSpline` to the **Rcpp** interface:
  [issue 19](https://github.com/wenjie2wang/splines2/issues/19).
* Added a new argument named `coef` to the `predict()` methods to compute the
  responding spline function and made it possible to obtain the derivatives or
  update spline basis functions by passing `...` to the `update()` methods.
* Added a new argument named `trim` to `naturalSpline()` to set the default
  boundary knots after trimming a fraction of observations.
* Added a new argument named `warn.outside` and a package option named
  `splines2.warn.outside` to specify if a warning should be thrown out for
  B-splines, etc. when any `x` is placed outside the boundary.
* Added the following function aliases to encourage the use in a model formula:
  * `bsp()` = `bSpline()`
  * `msp()` = `mSpline()`
  * `isp()` = `iSpline()`
  * `csp()` = `cSpline()`
  * `nsp()` = `naturalSpline()`
  * `bpoly()` = `bernsteinPoly()`
* Added a matrix named `H` to the attribution of objects for natural cubic
  splines so that users may transform cubic B-splines (from other
  software/packages) to the natural cubic splines (returned by
  `naturalSpline()`/`nsp()` or `nsk()`).
* Added an integer vector named `x_index` to the attribution of objects for
  regression splines (except Bernstein polynomials).


## Major changes

* Adjusted the class order of the returned objects.
* Adjusted the default placement of the internal knots from the specified `df`
  to be equidistant if the internal knots resulted from quantiles are
  problematic.  A warning will be thrown out in that case.


# splines2 0.4.8

## Bug fixes

* Fixed the Rcpp interface of `PeriodicMSpline` so that a simple knot sequence
  can be specified through `set_knot_sequence`:
  [issue 18](https://github.com/wenjie2wang/splines2/issues/18).


# splines2 0.4.7

## Minor changes

* Adjusted the column arrangement of the natural cubic spline basis matrix so
  that it matches the equations given in the JDS paper:
  [issue 17](https://github.com/wenjie2wang/splines2/issues/17).


# splines2 0.4.6

## New features

* Added `update()` methods to produce new spline basis functions based on the
  given object with specified updates in terms of `degree` and `knots`, etc.

## Minor changes

* Appended a new class named `splines2` to the output matrices to simplify some
  common S3 methods.


# splines2 0.4.5

## Minor changes

* Improved the computational efficiency of finding the knot intervals for `x`
  (by replacing the naive binary search implementation with `std::upper_bound`
  and `std::distance`).


# splines2 0.4.4

## New features

* Added the `makepredictcall()` methods for all available spline basis functions
  to help `model.frame.default()` create the right matrices when predicting from
  models with terms such as `bSpline()`, etc.  Thanks Zheyuan Li for suggesting
  this feature.
* Added arguments `derivs` and `integal` to `bSpline()` for consistency with
  `mSpline()` and `bernsteinPoly()`, etc.

## Minor changes

* Made the internal checking procedure more strict to throw an error if any
  internal knots are placed at or outside boundary:
  [issue 5](https://github.com/wenjie2wang/splines2/issues/5).

## Bug fixes

* Fixed the `predict()` method for `cSpline` objects when `scale = FALSE`.


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

* Updated the generation of the knot sequence for periodic M-splines following
  Piegl and Tiller (1997), which relaxed the previous requirement that
  `length(knots) >= degree` to `length(knots) >= degree - 1`.


# splines2 0.4.1

## New features

* Added function `naturalSpline()` providing implementation of nonnegative
  natural cubic splines.
* Added argument `periodic` to function `mSpline()` for periodic M-splines.
* Added argument `integral` to function `mSpline()` for integrals of M-splines
  or periodic M-splines.
* Added `deriv()`, `predict()`, and `print()` method for `naturalSpline` class
  object.

## Minor changes

* Updated the `deriv()` method for `mSpline` class object for periodic
  M-splines.


# splines2 0.3.1

## Minor changes

* Modified testing examples for CRAN tests on r-patched-solaris-x86 and
  r-oldrel-macos-x86_64.


# splines2 0.3.0

## New features

* Added function `bernsteinPoly()` providing implementation of generalized
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
