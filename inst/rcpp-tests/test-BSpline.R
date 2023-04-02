Rcpp::sourceCpp("test-BSpline.cpp")

x <- seq.int(0, 10, 0.02)
inter_knots <- c(2.4, 3.5, 5.2, 8)
bound_knots <- c(- 1, 12)
degree <- 4

foo <- function(...) {
    mat <- bSpline(..., intercept = TRUE)
    imat <- ibs(..., intercept = TRUE)
    d1mat <- deriv(mat)
    d2mat <- deriv(d1mat)
    d3mat <- deriv(d2mat)
    list(basis = unclass(mat),
         integral = unclass(imat),
         d1 = unclass(d1mat),
         d2 = unclass(d2mat),
         d3 = unclass(d3mat),
         degree = attr(mat, "degree"),
         internal_knots = knots(mat),
         boundary_knots = attr(mat, "Boundary.knots"))
}

res <- foo(x = x, knots = inter_knots, degree = degree,
           Boundary.knots = bound_knots)

## default constructors with setter methods
res00 <- rcpp_bspline00(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res00)

res01 <- rcpp_bspline01(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res01)

res02 <- rcpp_bspline02(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res02)

res03 <- rcpp_bspline03(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res03)

res04 <- rcpp_bspline04(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res04)

res05 <- rcpp_bspline05(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res05)

## non-default constructor 1
res1 <- rcpp_bspline1(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res1)

## non-default constructor 2
res2 <- rcpp_bspline2(x, 10, degree, bound_knots)
res20 <- foo(x = x, degree = degree, df = 10,
             Boundary.knots = bound_knots)
expect_equivalent(res20, res2)

## non-default constructor 3: simple knot sequence
knot_seq <- sort(c(rep(bound_knots, each = degree + 1), inter_knots))
res31 <- rcpp_bspline3(x, degree, knot_seq)
expect_equivalent(res, res31)

## non-default constructor 3: extended knot sequence
knot_seq <- sort(c(seq.int(0, 10, 1), 1, rep(4, 3), rep(7, 2)))
res32 <- rcpp_bspline3(x, degree, knot_seq)
expect_equivalent(
    res32$basis,
    splines::splineDesign(knot_seq, x, ord = degree + 1, outer.ok = TRUE)
)
expect_equivalent(
    res32$d1,
    splines::splineDesign(knot_seq, x, ord = degree + 1, outer.ok = TRUE,
                          derivs = 1)
)
expect_equivalent(
    res32$d2,
    splines::splineDesign(knot_seq, x, ord = degree + 1, outer.ok = TRUE,
                          derivs = 2)
)

## non-default constructor 4
res4 <- rcpp_bspline4(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res4)

## conversion from BernsteinPoly
res5 <- rcpp_bspline5(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res5)

## conversion from PeriodicMSpline
res6 <- rcpp_bspline6(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res6)
