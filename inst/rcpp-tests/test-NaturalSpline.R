Rcpp::sourceCpp("test-NaturalSpline.cpp")

x <- seq.int(0, 10, 0.02)
inter_knots <- c(2.4, 3.5, 5.2, 8)
bound_knots <- c(- 1, 12)

foo <- function(...) {
    mat <- naturalSpline(..., intercept = TRUE)
    imat <- naturalSpline(..., intercept = TRUE, integral = TRUE)
    d1mat <- deriv(mat)
    d2mat <- deriv(d1mat)
    d3mat <- deriv(d2mat)
    list(basis = mat, integral = imat,
         d1 = d1mat, d2 = d2mat, d3 = d3mat,
         degree = 3,
         internal_knots = knots(mat),
         boundary_knots = attr(mat, "Boundary.knots"))
}

res <- foo(x = x, knots = inter_knots, Boundary.knots = bound_knots)

## default constructors with setter methods
res00 <- rcpp_nspline00(x, inter_knots, bound_knots)
expect_equivalent(res, res00)

res01 <- rcpp_nspline01(x, inter_knots, bound_knots)
expect_equivalent(res, res01)

res02 <- rcpp_nspline02(x, inter_knots, bound_knots)
expect_equivalent(res, res02)

res03 <- rcpp_nspline03(x, inter_knots, bound_knots)
expect_equivalent(res, res03)

res04 <- rcpp_nspline04(x, inter_knots, bound_knots)
expect_equivalent(res, res04)

res05 <- rcpp_nspline05(x, inter_knots, bound_knots)
expect_equivalent(res, res05)

## non-default constructor 1
res1 <- rcpp_nspline1(x, inter_knots, bound_knots)
expect_equivalent(res, res1)

## non-default constructor 2
res2 <- rcpp_nspline2(x, 10, bound_knots)
res20 <- foo(x = x, df = 10, Boundary.knots = bound_knots)
expect_equivalent(res20, res2)

## non-default constructor 4
res4 <- rcpp_nspline4(x, inter_knots, bound_knots)
expect_equivalent(res, res4)

## conversion from BernsteinPoly
res5 <- rcpp_nspline5(x, inter_knots, degree = 4, bound_knots)
expect_equivalent(res, res5)

## conversion from PeriodicNspline
res6 <- rcpp_nspline6(x, inter_knots, degree = 5, bound_knots)
expect_equivalent(res, res6)
