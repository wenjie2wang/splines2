## get implementations of v0.2.8 for reference
v2 <- new.env()
source("../v0.2.8.R", v2)

## helper functions
isNumMatrix <- v2$isNumMatrix

### 1. check correctness first
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
knots <- c(1, 3)
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

## compare with v0.2.8
x <- c(NA, seq.int(0, 0.5, 0.1), NA, seq.int(0.6, 1, 0.1), NA)
knots <- c(0.25, 0.5, 0.75)
x2 <- c(- 1, 2, x)
b_knots <- c(0, 1)

## default cubic splines without internal knots
expect_equivalent(ibs(x), v2$ibs(x))

## cubic splines with specified df
expect_equivalent(ibs(x, df = 5),
                  v2$ibs(x, df = 5))

## cubic splines with specified internal knots
expect_equivalent(ibs(x, knots = knots),
                  v2$ibs(x, knots = knots))

## qudractic splines without internal knots
expect_equivalent(ibs(x, degree = 2L),
                  v2$ibs(x, degree = 2L))

## complete basis with intercept
expect_equivalent(ibs(x, intercept = TRUE),
                  v2$ibs(x, intercept = TRUE))

## specified knots
expect_equivalent(ibs(x, knots = knots, intercept = TRUE),
                  v2$ibs(x, knots = knots, intercept = TRUE))

## specified df
expect_equivalent(ibs(x, df = 6, intercept = TRUE),
                  v2$ibs(x, df = 6, intercept = TRUE))

## degree zero
expect_equivalent(ibs(x, df = 5, degree = 0),
                  v2$ibs(x, df = 5, degree = 0))
expect_equivalent(ibs(x, df = 5, degree = 0, intercept = TRUE),
                  v2$ibs(x, df = 5, degree = 0, intercept = TRUE))
bsMat0a <- ibs(x, degree = 0, intercept = TRUE)
bsMat0b <- ibs(x, df = 5, degree = 0)
bsMat0c <- ibs(x, df = 5, degree = 0, intercept = TRUE)
bsMat0d <- ibs(x, knots = knots, degree = 0)
bsMat0e <- ibs(x, knots = knots, degree = 0, intercept = TRUE)
expect_true(isNumMatrix(bsMat0a, 14L, 1L))
expect_equal(sum(is.na(bsMat0b)), 15L) # keep NA's as is
expect_true(isNumMatrix(bsMat0b, 14L, 5L))
expect_true(isNumMatrix(bsMat0c, 14L, 5L))
expect_true(isNumMatrix(bsMat0d, 14L, 3L))
expect_true(isNumMatrix(bsMat0e, 14L, 4L))
expect_true(isNumMatrix(
    ibs(x, df = 10, knots = knots, degree = 0L),
    14L, 3L))
expect_true(isNumMatrix(
    ibs(x, df = 10, knots = knots,
        degree = 0, intercept = TRUE),
    14L, 4L))

## x outside of boundary
suppressWarnings({
    expect_equivalent(
        ibs(x2, df = 6, degree = 3, Boundary.knots = b_knots),
        v2$ibs(x2, df = 6, degree = 3, Boundary.knots = b_knots)
    )
})
suppressWarnings({
    expect_equivalent(
        ibs(x2, knots = knots, degree = 3, Boundary.knots = b_knots),
        v2$ibs(x2, knots = knots, degree = 3, Boundary.knots = b_knots)
    )
})

## keep names of x
names(x) <- sample(LETTERS, length(x), replace = TRUE)
expect_equal(rownames(ibs(x)), names(x))

### 2. check designed features with expectation
## NA is only allowed in x

## error if all of x are NA's
expect_error(ibs(c(NA_real_, NA_real_), degree = 0))
expect_error(ibs(c(NA, NA), df = 5))

## error if degree has NA or negative
expect_error(ibs(x, degree = - 1))
expect_error(ibs(x, degree = c(2, NA)))
expect_error(ibs(x, degree = NA))

## error if df has NA or negative
expect_error(ibs(x, df = - 1))
expect_error(ibs(x, df = c(2, NA)))
expect_error(ibs(x, df = NA))

## error if knots has NA
expect_error(ibs(x, knots = c(0.1, 0.5, NA)), "knots")
expect_error(ibs(x, Boundary.knots = c(0.1, 0.5, NA)), "knots")

## error if boundary knots are inappropriate
expect_error(ibs(x, Boundary.knots = 0.1), "knots")
expect_error(ibs(x, Boundary.knots = c(0.1, 0.1)), "knots")
expect_error(ibs(x, Boundary.knots = c(0.1, 0.5, 1)), "knots")

## error if empty matrix
expect_true(isNumMatrix(ibs(x, degree = 0, intercept = TRUE),
                        length(x), 1))
expect_error(ibs(x, degree = 0), "No column")

## error if any internal knot is not placed inside of boundary
expect_error(ibs(x, knots = c(0, 0.5), degree = 0), "inside")

## warning if any x outside of boundary
expect_warning(ibs(c(x, 10), knots = knots, degree = 0,
                   Boundary.knots = c(0, 1)),
               "beyond boundary knots")
expect_warning(ibs(c(x, 10), knots = knots, degree = 3,
                   Boundary.knots = c(0, 1)),
               "beyond boundary knots")
