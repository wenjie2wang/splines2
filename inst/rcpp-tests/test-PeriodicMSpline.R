Rcpp::sourceCpp("test-PeriodicMSpline.cpp")

x <- seq.int(0, 10, 0.02)
inter_knots <- c(2.4, 3.5, 5.2, 8)
bound_knots <- c(- 1, 12)
degree <- 4
order <- degree + 1
knot_seq <- c(rep(bound_knots[1], order), inter_knots,
              rep(bound_knots[2], order))

foo <- function(...) {
    mat <- mSpline(..., intercept = TRUE, periodic = TRUE)
    imat <- mSpline(..., intercept = TRUE, periodic = TRUE, integral = TRUE)
    d1mat <- deriv(mat)
    d2mat <- deriv(d1mat)
    d3mat <- deriv(d2mat)
    list(basis = mat, integral = imat,
         d1 = d1mat, d2 = d2mat, d3 = d3mat,
         degree = attr(mat, "degree"),
         internal_knots = knots(mat),
         boundary_knots = attr(mat, "Boundary.knots"))
}

res <- foo(x = x, knots = inter_knots, degree = degree,
           Boundary.knots = bound_knots)

## default constructors with setter methods
res00 <- rcpp_pmspline00(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res00)

res01 <- rcpp_pmspline01(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res01)

res02 <- rcpp_pmspline02(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res02)

res03 <- rcpp_pmspline03(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res03)

res04 <- rcpp_pmspline04(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res04)

res05 <- rcpp_pmspline05(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res05)

res06 <- rcpp_pmspline06(x, degree, knot_seq)
expect_equivalent(res, res06)

res07 <- rcpp_pmspline07(x, degree, knot_seq)
expect_equivalent(res, res07)

res08 <- rcpp_pmspline08(x, degree, knot_seq)
expect_equivalent(res, res08)

## non-default constructor 1
res1 <- rcpp_pmspline1(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res1)

## non-default constructor 2
res2 <- rcpp_pmspline2(x, 10, degree, bound_knots)
res20 <- foo(x = x, degree = degree, df = 10,
             Boundary.knots = bound_knots)
expect_equivalent(res20, res2)

## non-default constructor 3
res3 <- rcpp_pmspline3(x, degree, knot_seq)
expect_equivalent(res, res3)

## non-default constructor 4
res4 <- rcpp_pmspline4(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res4)

## conversion from BernsteinPoly
res5 <- rcpp_pmspline5(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res5)
