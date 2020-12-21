## get implementations of v0.2.8 for reference
v2 <- new.env()
source("../v0.2.8.R", v2)

## helper functions
isNumMatrix <- v2$isNumMatrix

### 1. check correctness first
x <- seq.int(0, 1, 0.1)
knots <- c(0.3, 0.5, 0.6)
b_knots <- c(0, 1)

## natural spline basis
nsMat0 <- naturalSpline(x, knots = knots, intercept = TRUE)
## integrals
nsMat1 <- naturalSpline(x, knots = knots, intercept = TRUE, integral = TRUE)
## first derivatives
nsMat2 <- naturalSpline(x, knots = knots, intercept = TRUE, derivs = 1)
## second derivatives
nsMat3 <- naturalSpline(x, knots = knots, intercept = TRUE, derivs = 2)

## check matrix size
expect_true(isNumMatrix(nsMat0, length(x), length(knots) + 2L))
expect_true(isNumMatrix(nsMat1, length(x), length(knots) + 2L))
expect_true(isNumMatrix(nsMat2, length(x), length(knots) + 2L))
expect_true(isNumMatrix(nsMat3, length(x), length(knots) + 2L))

## specify df directly instead of knots
for (j in seq.int(2, 10)) {
    expect_true(isNumMatrix(
        naturalSpline(x, df = j), length(x), j
    ))
}

## use the deriv method
expect_equivalent(nsMat0, deriv(nsMat1))
expect_equivalent(nsMat2, deriv(nsMat0))
expect_equivalent(nsMat3, deriv(nsMat2))
expect_equivalent(nsMat3, deriv(nsMat0, 2))

## check second derivatives at boundary knots
expect_true(all(abs(predict(nsMat3, 0)) < 1e-12))
expect_true(all(abs(predict(nsMat3, 1)) < 1e-12))

### 2. checking inputs
x <- c(NA, seq.int(0, 0.5, 0.1), NA, seq.int(0.6, 1, 0.1), NA)
knots <- c(0.25, 0.5, 0.8)
x2 <- c(- 1, 2, x)
b_knots <- c(0, 1)

## warning for having x outside of boundary
expect_warning(naturalSpline(x2, df = 5, Boundary.knots = b_knots))
