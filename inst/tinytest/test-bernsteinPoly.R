source("utils.R")

### 1. test correctness first
x <- seq.int(0, 1, 0.1)

## degree 0: basis
res0 <- bernsteinPoly(x, degree = 0, intercept = TRUE)
expect_eqt(matrix(1, nrow = length(x)), res0)
res0 <- bernsteinPoly(x, degree = 0, intercept = TRUE,
                      derivs = 1, integral = TRUE)
expect_eqt(matrix(1, nrow = length(x)), res0)
## degree 0: derivative
res0 <- bernsteinPoly(x, degree = 0, intercept = TRUE, derivs = 1)
expect_eqt(matrix(0, nrow = length(x)), res0)
## degree 0: integral
res0 <- bernsteinPoly(x, degree = 0, intercept = TRUE, integral = TRUE)
expect_eqt(matrix(x, nrow = length(x)), res0)


## degree 1: basis
bp1 <- function(x) { cbind(1 - x, x) }
res1 <- bernsteinPoly(x, degree = 1, intercept = TRUE)
expect_eqt(bp1(x), res1)
res1 <- bernsteinPoly(x, degree = 1, intercept = FALSE,
                      derivs = 1, integral = TRUE)
expect_eqt(bp1(x)[, -1L, drop = FALSE], res1)
## degree 1: first derivative
dbp1 <- function(x) {
    matrix(c(- 1, 1), nrow = length(x), ncol = 2, byrow = TRUE)
}
res1 <- bernsteinPoly(x, degree = 1, intercept = TRUE, derivs = 1)
expect_eqt(dbp1(x), res1)
res1 <- bernsteinPoly(x, degree = 1, intercept = FALSE,
                      derivs = 2, integral = TRUE)
expect_eqt(dbp1(x)[, - 1L, drop = FALSE], res1)
## degree 1: second derivative
ddbp1 <- function(x) { matrix(0, nrow = length(x), ncol = 2) }
res1 <- bernsteinPoly(x, degree = 1, intercept = TRUE, derivs = 2)
expect_eqt(ddbp1(x), res1)
res1 <- bernsteinPoly(x, degree = 1, intercept = FALSE,
                      derivs = 3, integral = TRUE)
expect_eqt(ddbp1(x)[, - 1L, drop = FALSE], res1)
## degree 1: integral
ibp1 <- function(x) {
    cbind(x - 0.5 * x ^ 2, 0.5 * x ^ 2)
}
res1 <- bernsteinPoly(x, degree = 1, intercept = TRUE, integral = TRUE)
expect_eqt(ibp1(x), res1)
res1 <- bernsteinPoly(x, degree = 1, intercept = FALSE, integral = TRUE)
expect_eqt(ibp1(x)[, - 1L, drop = FALSE], res1)

## degree 2: basis
bp2 <- function(x) {
    cbind((1 - x) ^ 2, 2 * x * (1 - x), x ^ 2)
}
res2 <- bernsteinPoly(x, degree = 2, intercept = TRUE)
expect_eqt(bp2(x), res2)
## degree 2: first derivative
dbp2 <- function(x) {
    cbind(- 2 * (1 - x), 2 * (1 - 2 * x), 2 * x)
}
res2 <- bernsteinPoly(x, degree = 2, intercept = TRUE, derivs = 1)
expect_eqt(dbp2(x), res2)
res2 <- bernsteinPoly(x, degree = 2, intercept = TRUE,
                      derivs = 2, integral = TRUE)
expect_eqt(dbp2(x), res2)
## degree 2: second derivative
ddbp2 <- function(x) {
    matrix(c(2, - 4, 2), nrow = length(x), ncol = 3, byrow = TRUE)
}
res2 <- bernsteinPoly(x, degree = 2, intercept = TRUE,
                      derivs = 2, integral = FALSE)
expect_eqt(ddbp2(x), res2)
## degree 2: integral
ibp2 <- function(x) {
    cbind(1 / 3 - (1 - x) ^ 3 / 3, x ^ 2 - 2 / 3 * x ^ 3, 1 / 3 * x ^ 3)
}
res2 <- bernsteinPoly(x, degree = 2, intercept = TRUE, integral = TRUE)
expect_eqt(ibp2(x), res2)

## degree 3: basis
bp3 <- function(x) {
    cbind((1 - x) ^ 3, 3 * x * (1 - x) ^ 2, 3 * (1 - x) * x ^ 2, x ^ 3)
}
res3 <- bernsteinPoly(x, degree = 3, intercept = TRUE)
expect_eqt(bp3(x), res3)
## degree 3: first derivative
dbp3 <- function(x) {
    cbind(- 3 * (1 - x) ^ 2, 3 * ((1 - x) ^ 2 - 2 * x * (1 - x)),
          3 * (- x ^ 2 + (1 - x) * 2 * x), 3 * x ^ 2)
}
res3 <- bernsteinPoly(x, degree = 3, intercept = TRUE, derivs = 1)
expect_eqt(dbp3(x), res3)
res3 <- bernsteinPoly(x, degree = 3, intercept = TRUE,
                      derivs = 2, integral = TRUE)
expect_eqt(dbp3(x), res3)
## degree 3: second derivative
ddbp3 <- function(x) {
    cbind(6 * (1 - x), 3 * (6 * x - 4), 3 * (- 6 * x + 2), 6 * x)
}
res3 <- bernsteinPoly(x, degree = 3, intercept = TRUE, derivs = 2)
expect_eqt(ddbp3(x), res3)
## degree 3: integral
ibp3 <- function(x) {
    cbind(1 / 4 - (1 - x) ^ 4 / 4,
          3 / 2 * x ^ 2 - 2 * x ^ 3 + 3 / 4 * x ^ 4,
          x ^ 3 - 3 / 4 * x ^ 4, x ^ 4 / 4)
}
res3 <- bernsteinPoly(x, degree = 3, intercept = TRUE, integral = TRUE)
expect_eqt(ibp3(x), res3)

## for given boundary knots
a <- - 1.23
b <- 2.68
get_x <- function(x, a, b) {
    (x - a) / (b - a)
}
xx <- seq.int(a, b, by = 0.01)

## basis
res0 <- bp3(get_x(xx, a, b))
res1 <- bernsteinPoly(xx, degree = 3, intercept = TRUE,
                      Boundary.knots = c(a, b))
expect_eqt(res0, res1)

## first derivatives
res0 <- dbp3(get_x(xx, a, b)) / (b - a)
res1 <- bernsteinPoly(xx, degree = 3, intercept = TRUE,
                      Boundary.knots = c(a, b), derivs = 1)
expect_eqt(res0, res1)

## second derivatives
res0 <- ddbp3(get_x(xx, a, b)) / (b - a) ^ 2
res1 <- bernsteinPoly(xx, degree = 3, intercept = TRUE,
                      Boundary.knots = c(a, b), derivs = 2)
expect_eqt(res0, res1)

## integral
res0 <- ibp3(get_x(xx, a, b)) * (b - a)
res1 <- bernsteinPoly(xx, degree = 3, intercept = TRUE,
                      Boundary.knots = c(a, b), integral = TRUE)
expect_eqt(res0, res1)

## for named x
named_x <- seq.int(- 1, 1, 0.1)
names(named_x) <- paste("x =", sprintf("%.1f", named_x))
res3 <- bernsteinPoly(named_x, degree = 5)
expect_equal(row.names(res3), names(named_x))

### 2. edge cases
## NA is only allowed in x
## error if all of x are NA's
expect_error(bernsteinPoly(c(NA_real_, NA_real_), degree = 0))
expect_error(bernsteinPoly(c(NA, NA), degree = 5))

## error if degree has NA or negative
expect_error(bernsteinPoly(x, degree = - 1))
expect_error(bernsteinPoly(x, degree = NA))

## error if empty matrix
expect_error(bernsteinPoly(x, degree = 0))

## error if x is outside of boundary
expect_error(bernsteinPoly(c(- 1, x, 1, 1.1), degree = 2,
                           Boundary.knots = c(0, 1)))

## error if the specified boundary knots is inappropriate
expect_error(bernsteinPoly(c(- 1, x, 1, 1.1), degree = 2,
                           Boundary.knots = 0))
