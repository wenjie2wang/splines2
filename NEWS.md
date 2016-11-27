# CHANGES IN splines2 VERSION 0.2.2.9000

## BUG FIXES

* Fixed one-piece constant basis for M-splines.


# CHANGES IN splines2 VERSION 0.2.2

## BUG FIXES

* Fixed the NA's handling in all the functions constructing spline bases.


# CHANGES IN splines2 VERSION 0.2.1

## NEW FEATURES

* Added function `bSpline` generating B-spline basis allowing zero degree or
  piecewise constant basis based on function `bs` in package **splines**.

* Introduced function `bSpline` to allow M-splines of degree zero.

* Added function `cSpline` constructing convex spline (C-spline) basis.

* Added `predict` methods for `bSpline2` object and `cSpline` object generated
  by `bSpline` and `cSpline`, respectively.

* Added `print` methods for all **splines2** objects developed so far.

## MAJOR CHANGES

* Improved the function `iSpline` to construct I-spline basis directly from
  B-spline basis instead of M-spline basis.

## MINOR CHANGES

* Updated all CRAN URL to a canonical form suggested.


# CHANGES IN splines2 VERSION 0.1.0

## NEW FEATURES

* The first version of **splines2** providing functions constructing M-spline,
  I-spline, and integral of B-spline basis.


