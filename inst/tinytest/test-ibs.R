## without internal knots
x <- seq.int(0, 1, 0.1)
## degree = 0
expect_equivalent(matrix(x), ibs(x, degree = 0, intercept = TRUE))
## degree = 1
b1 <- function(x) x - x ^ 2 / 2
b2 <- function(x) x ^ 2 / 2
i1mat <- cbind(b1(x), b2(x))
expect_equivalent(matrix(b2(x)), ibs(x, degree = 1))
expect_equivalent(i1mat, ibs(x, degree = 1, intercept = TRUE))
## degree = 2
b1 <- function(x) x ^ 3 / 3 - x ^ 2 + x
b2 <- function(x) - 2 * x ^ 3 / 3 + x ^ 2
b3 <- function(x) x ^ 3 / 3
i2mat <- cbind(b1(x), b2(x), b3(x))
expect_equivalent(cbind(b2(x), b3(x)), ibs(x, degree = 2))
expect_equivalent(i2mat, ibs(x, degree = 2, intercept = TRUE))
## degree = 3
b1 <- function(x) - (1 - x) ^ 4 / 4 + 1 / 4
b2 <- function(x) 3 / 4 * x ^ 4 - 2 * x ^ 3 + 3 / 2 * x ^ 2
b3 <- function(x) x ^ 3 - 3 / 4 * x ^ 4
b4 <- function(x) x ^ 4 / 4
i3mat <- cbind(b1(x), b2(x), b3(x), b4(x))
expect_equivalent(cbind(b2(x), b3(x), b4(x)), ibs(x, degree = 3))
expect_equivalent(i3mat, ibs(x, degree = 3, intercept = TRUE))

## with two internal knots
x <- seq.int(0, 4, 0.1)
knots <- c(NA, 1, NA, 3, NA)
ind01 <- function(x) as.numeric(x >= 0 & x < 1)
ind13 <- function(x) as.numeric(x >= 1 & x < 3)
ind34 <- function(x) as.numeric(x >= 3 & x <= 4)
## degree = 0
b1 <- function(x) ind01(x) * x + ind13(x) + ind34(x)
b2 <- function(x) ind13(x) * (x - 1) + 2 * ind34(x)
b3 <- function(x) ind34(x) * (x - 3)
expect_equivalent(cbind(b2(x), b3(x)),
                  ibs(x, knots = knots, degree = 0))
expect_equivalent(cbind(b1(x), b2(x), b3(x)),
                  ibs(x, knots = knots, degree = 0, intercept = TRUE))
## degree = 1
b1 <- function(x) ind01(x) * (x - x ^ 2 / 2) + (ind13(x) + ind34(x)) / 2
b2 <- function(x) {
    ind01(x) * x ^ 2 / 2 +
        ind13(x) * (1.5 * x - 0.25 * x ^ 2 - 3 / 4) + ind34(x) * 1.5
}
b3 <- function(x) {
    ind13(x) * (x ^ 2 / 4 - x / 2 + 1 / 4) +
        ind34(x) * (- x ^ 2 / 2 + 4 * x - 6.5)
}
b4 <- function(x) ind34(x) * (x ^ 2 / 2 - 3 * x + 4.5)
expect_equivalent(cbind(b2(x), b3(x), b4(x)),
                  ibs(x, knots = knots, degree = 1))
expect_equivalent(cbind(b1(x), b2(x), b3(x), b4(x)),
                  ibs(x, knots = knots, degree = 1, intercept = TRUE))
