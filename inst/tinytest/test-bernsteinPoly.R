### 1. test correctness first
x <- seq.int(0, 1, 0.1)

## degree 0: basis
res0 <- bernsteinPoly(x, degree = 0, intercept = TRUE)
expect_equivalent(matrix(1, nrow = length(x)), res0)
res0 <- bernsteinPoly(x, degree = 0, intercept = TRUE,
                      derivs = 1, integral = TRUE)
expect_equivalent(matrix(1, nrow = length(x)), res0)
## degree 0: derivative
res0 <- bernsteinPoly(x, degree = 0, intercept = TRUE, derivs = 1)
expect_equivalent(matrix(0, nrow = length(x)), res0)
## degree 0: integral
res0 <- bernsteinPoly(x, degree = 0, intercept = TRUE, integral = TRUE)
expect_equivalent(matrix(x, nrow = length(x)), res0)


## degree 1: basis
bp1 <- function(x) { cbind(1 - x, x) }
res1 <- bernsteinPoly(x, degree = 1, intercept = TRUE)
expect_equivalent(bp1(x), res1)
res1 <- bernsteinPoly(x, degree = 1, intercept = FALSE,
                      derivs = 1, integral = TRUE)
expect_equivalent(bp1(x)[, -1L, drop = FALSE], res1)
## degree 1: first derivative
dbp1 <- function(x) {
    matrix(c(- 1, 1), nrow = length(x), ncol = 2, byrow = TRUE)
}
res1 <- bernsteinPoly(x, degree = 1, intercept = TRUE, derivs = 1)
expect_equivalent(dbp1(x), res1)
res1 <- bernsteinPoly(x, degree = 1, intercept = FALSE,
                      derivs = 2, integral = TRUE)
expect_equivalent(dbp1(x)[, - 1L, drop = FALSE], res1)
## degree 1: second derivative
ddbp1 <- function(x) { matrix(0, nrow = length(x), ncol = 2) }
res1 <- bernsteinPoly(x, degree = 1, intercept = TRUE, derivs = 2)
expect_equivalent(ddbp1(x), res1)
res1 <- bernsteinPoly(x, degree = 1, intercept = FALSE,
                      derivs = 3, integral = TRUE)
expect_equivalent(ddbp1(x)[, - 1L, drop = FALSE], res1)
## degree 1: integral
ibp1 <- function(x) {
    cbind(x - 0.5 * x ^ 2, 0.5 * x ^ 2)
}
res1 <- bernsteinPoly(x, degree = 1, intercept = TRUE, integral = TRUE)
expect_equivalent(ibp1(x), res1)
res1 <- bernsteinPoly(x, degree = 1, intercept = FALSE, integral = TRUE)
expect_equivalent(ibp1(x)[, - 1L, drop = FALSE], res1)

## degree 2: basis
bp2 <- function(x) {
    cbind((1 - x) ^ 2, 2 * x * (1 - x), x ^ 2)
}
res2 <- bernsteinPoly(x, degree = 2, intercept = TRUE)
expect_equivalent(bp2(x), res2)
## degree 2: first derivative
dbp2 <- function(x) {
    cbind(- 2 * (1 - x), 2 * (1 - 2 * x), 2 * x)
}
res2 <- bernsteinPoly(x, degree = 2, intercept = TRUE, derivs = 1)
expect_equivalent(dbp2(x), res2)
res2 <- bernsteinPoly(x, degree = 2, intercept = TRUE,
                      derivs = 2, integral = TRUE)
expect_equivalent(dbp2(x), res2)
## degree 2: second derivative
ddbp2 <- function(x) {
    matrix(c(2, - 4, 2), nrow = length(x), ncol = 3, byrow = TRUE)
}
res2 <- bernsteinPoly(x, degree = 2, intercept = TRUE,
                      derivs = 2, integral = FALSE)
expect_equivalent(ddbp2(x), res2)
## degree 2: integral
ibp2 <- function(x) {
    cbind(1 / 3 - (1 - x) ^ 3 / 3, x ^ 2 - 2 / 3 * x ^ 3, 1 / 3 * x ^ 3)
}
res2 <- bernsteinPoly(x, degree = 2, intercept = TRUE, integral = TRUE)
expect_equivalent(ibp2(x), res2)

## degree 3: basis
bp3 <- function(x) {
    cbind((1 - x) ^ 3, 3 * x * (1 - x) ^ 2, 3 * (1 - x) * x ^ 2, x ^ 3)
}
res3 <- bernsteinPoly(x, degree = 3, intercept = TRUE)
expect_equivalent(bp3(x), res3)
## degree 3: first derivative
dbp3 <- function(x) {
    cbind(- 3 * (1 - x) ^ 2, 3 * ((1 - x) ^ 2 - 2 * x * (1 - x)),
          3 * (- x ^ 2 + (1 - x) * 2 * x), 3 * x ^ 2)
}
res3 <- bernsteinPoly(x, degree = 3, intercept = TRUE, derivs = 1)
expect_equivalent(dbp3(x), res3)
res3 <- bernsteinPoly(x, degree = 3, intercept = TRUE,
                      derivs = 2, integral = TRUE)
expect_equivalent(dbp3(x), res3)
## degree 3: second derivative
ddbp3 <- function(x) {
    cbind(6 * (1 - x), 3 * (6 * x - 4), 3 * (- 6 * x + 2), 6 * x)
}
res3 <- bernsteinPoly(x, degree = 3, intercept = TRUE, derivs = 2)
expect_equivalent(ddbp3(x), res3)
## degree 3: integral
ibp3 <- function(x) {
    cbind(1 / 4 - (1 - x) ^ 4 / 4,
          3 / 2 * x ^ 2 - 2 * x ^ 3 + 3 / 4 * x ^ 4,
          x ^ 3 - 3 / 4 * x ^ 4, x ^ 4 / 4)
}
res3 <- bernsteinPoly(x, degree = 3, intercept = TRUE, integral = TRUE)
expect_equivalent(ibp3(x), res3)


### 2. edge cases
## NA is only allowed in x
## error if all of x are NA's
expect_error(bernsteinPoly(c(NA_real_, NA_real_), degree = 0))
expect_error(bernsteinPoly(c(NA, NA), degree = 5))

## error if degree has NA or negative
expect_error(bernsteinPoly(x, degree = - 1))
expect_error(bernsteinPoly(x, degree = NA))

## error if empty matrix
expect_error(bernsteinPoly(x, degree = 0), "No column")

## error if x is outside of [0, 1]
expect_error(bernsteinPoly(c(- 1, x, 1, 1.1), degree = 2), "[0, 1]")
