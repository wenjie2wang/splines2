Rcpp::sourceCpp("test-BernsteinPoly.cpp")

x <- seq.int(0, 10, 0.02)
inter_knots <- c(2.4, 3.5, 5.2, 8)
bound_knots <- c(- 1, 12)
degree <- 4

foo <- function(...) {
    mat <- bernsteinPoly(..., intercept = TRUE)
    imat <- bernsteinPoly(..., intercept = TRUE, integral = TRUE)
    d1mat <- deriv(mat)
    d2mat <- deriv(d1mat)
    d3mat <- deriv(d2mat)
    list(basis = unclass(mat),
         integral = unclass(imat),
         d1 = unclass(d1mat),
         d2 = unclass(d2mat),
         d3 = unclass(d3mat),
         degree = attr(mat, "degree"),
         boundary_knots = attr(mat, "Boundary.knots"))
}

res <- foo(x = x, degree = degree,
           Boundary.knots = bound_knots)

## default constructors with setter methods
res00 <- rcpp_bp00(x, degree, bound_knots)
expect_equivalent(res, res00)

res01 <- rcpp_bp01(x, degree, bound_knots)
expect_equivalent(res, res01)

res02 <- rcpp_bp02(x, degree, bound_knots)
expect_equivalent(res, res02)

res03 <- rcpp_bp03(x, degree, bound_knots)
expect_equivalent(res, res03)

res04 <- rcpp_bp04(x, degree, bound_knots)
expect_equivalent(res, res04)

res05 <- rcpp_bp05(x, degree, bound_knots)
expect_equivalent(res, res05)

## non-default constructor 1
res1 <- rcpp_bp1(x, degree, bound_knots)
expect_equivalent(res, res1)

## non-default constructor 4
res4 <- rcpp_bp4(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res4)

## conversion from PeriodicMSpline
res6 <- rcpp_bp6(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res6)
