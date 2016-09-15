# splines2 v0.1.0.9000

## New features

* Added function `bSpline` generating B-spline basis allowing zero degree or
  piecewise constant basis based on function `bs` in package *splines*.

* Added function `cSpline` constructing convex spline (C-spline) basis.

* Added `predict` methods for `bSpline2` object and `cSpline` object generated
  by `bSpline` and `cSpline`, respectively.

## Major changes

* Improved the function `iSpline` that constructs I-spline basis directly from
  B-spline basis instead of M-spline basis.


# splines2 v0.1.0

## New features

* The first version of *splines2* providing functions constructing M-spline,
  I-spline, and integral of B-spline basis.


