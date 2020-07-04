## get implementations of v0.2.8 for reference
v2 <- new.env()
source("../v0.2.8.R", v2)

## helper functions
isNumMatrix <- v2$isNumMatrix

### 1. check correctness first
x <- c(NA, seq.int(0, 0.5, 0.1), NA, seq.int(0.6, 1, 0.1), NA)
knots <- c(0.25, 0.5, 0.75)
x2 <- c(- 1, 2, x)
b_knots <- c(0, 1)

## default cubic splines without internal knots
expect_equivalent(mSpline(x), v2$mSpline(x))

## cubic splines with specified df
expect_equivalent(mSpline(x, df = 5),
                  v2$mSpline(x, df = 5))

## cubic splines with specified internal knots
expect_equivalent(mSpline(x, knots = knots),
                  v2$mSpline(x, knots = knots))

## qudractic splines without internal knots
expect_equivalent(mSpline(x, degree = 2L),
                  v2$mSpline(x, degree = 2L))

## complete basis with intercept
expect_equivalent(mSpline(x, intercept = TRUE),
                  v2$mSpline(x, intercept = TRUE))

## specified knots
expect_equivalent(mSpline(x, knots = knots, intercept = TRUE),
                  v2$mSpline(x, knots = knots, intercept = TRUE))

## specified df
expect_equivalent(mSpline(x, df = 6, intercept = TRUE),
                  v2$mSpline(x, df = 6, intercept = TRUE))

## for testing splines with degree zero
msMat0a <- mSpline(x, degree = 0, intercept = TRUE)
msMat0b <- mSpline(x, knots = knots, degree = 1)
msMat0c <- mSpline(x, knots = knots, degree = 1, intercept = TRUE)
msMat0d <- mSpline(x, knots = knots, degree = 2)
msMat0e <- mSpline(x, knots = knots, degree = 2, intercept = TRUE)
msMat0f <- mSpline(0.1, knots = knots, degree = 2, intercept = TRUE,
                   Boundary.knots = c(0, 1), derivs = 1)
msMat0g <- mSpline(0.1, df = 6, degree = 2, intercept = TRUE,
                   Boundary.knots = c(0, 1), derivs = 2)
msMat0h <- mSpline(0.1, knots = knots, degree = 2,
                   Boundary.knots = c(0, 1), derivs = 3)
expect_true(isNumMatrix(msMat0a, 14L, 1L))
expect_equal(sum(is.na(msMat0b)), 12L) # keep NA's as is
expect_true(isNumMatrix(msMat0b, 14L, 4L))
expect_true(isNumMatrix(msMat0c, 14L, 5L))
expect_true(isNumMatrix(msMat0d, 14L, 5L))
expect_true(isNumMatrix(msMat0e, 14L, 6L))
expect_true(isNumMatrix(msMat0f, 1L, 6L))
expect_true(isNumMatrix(msMat0g, 1L, 6L))
expect_true(isNumMatrix(msMat0h, 1L, 5L))
expect_error(mSpline(x, degree = 0))
expect_warning(mSpline(c(x, 10), knots = knots, degree = 0,
                       Boundary.knots = c(0, 1)),
               "beyond boundary knots")

## true close form formula given the all knots and degree
## transformation of constant bases
x3 <- seq.int(0, 7, 0.1)
m0_1 <- function(x) as.numeric(x < 1)
m0_2 <- function(x) as.numeric(x >= 1 & x < 3) * 0.5
m0_3 <- function(x) as.numeric(x >= 3 & x <= 7) * 0.25
expect_equivalent(mSpline(x3, knots = c(1, 3), degree = 0L, intercept = TRUE),
                  cbind(m0_1(x3), m0_2(x3), m0_3(x3)))

## transformation of linear bases
x4 <- seq.int(0, 3, 0.1)
ind01 <- function(x) as.numeric(x >= 0 & x < 1)
ind12 <- function(x) as.numeric(x >= 1 & x < 2)
ind23 <- function(x) as.numeric(x >= 2 & x <= 3)
m0_1 <- function(x) ind01(x) * (2 - 2 * x)
m0_2 <- function(x) ind01(x) * x + ind12(x) * (2 - x)
m0_3 <- function(x) ind12(x) * (x - 1) + ind23(x) * (3 - x)
m0_4 <- function(x) ind23(x) * (2 * x - 4)
expect_equivalent(mSpline(x4, knots = c(1, 2), degree = 1L, intercept = TRUE),
                  cbind(m0_1(x4), m0_2(x4), m0_3(x4), m0_4(x4)))

## x outside of boundary
suppressWarnings({
    expect_equivalent(
        mSpline(x2, df = 6, degree = 3, Boundary.knots = b_knots),
        v2$mSpline(x2, df = 6, degree = 3, Boundary.knots = b_knots)
    )
})
suppressWarnings({
    expect_equivalent(
        mSpline(x2, knots = knots, degree = 3, Boundary.knots = b_knots),
        v2$mSpline(x2, knots = knots, degree = 3, Boundary.knots = b_knots)
    )
})

## keep names of x
names(x) <- sample(LETTERS, length(x), replace = TRUE)
expect_equal(rownames(bSpline(x)), names(x))

### 2. check designed features with expectation
## NA is only allowed in x

## error if all of x are NA's
expect_error(mSpline(c(NA_real_, NA_real_), degree = 0))
expect_error(mSpline(c(NA, NA), df = 5))

## error if degree has NA or negative
expect_error(mSpline(x, degree = - 1))
expect_error(mSpline(x, degree = NA))

## error if df has NA or negative
expect_error(mSpline(x, df = - 1))
expect_error(mSpline(x, df = NA))

## error if knots has NA
expect_error(mSpline(x, knots = c(0.1, 0.5, NA)), "knots")
expect_error(mSpline(x, Boundary.knots = c(0.1, 0.5, NA)), "knots")

## error if boundary knots are inappropriate
expect_error(mSpline(x, Boundary.knots = 0.1), "knots")
expect_error(mSpline(x, Boundary.knots = c(0.1, 0.1)), "knots")
expect_error(mSpline(x, Boundary.knots = c(0.1, 0.5, 1)), "knots")

## error if empty matrix
expect_true(isNumMatrix(mSpline(x, degree = 0, intercept = TRUE),
                        length(x), 1))
expect_error(mSpline(x, degree = 0), "No column")
expect_error(mSpline(x, degree = 0, derivs = 1), "No column")

## error if any internal knot is not placed inside of boundary
expect_error(mSpline(x, knots = c(0, 0.5), degree = 0), "inside")

## warning if any x outside of boundary
expect_warning(mSpline(c(x, 10), knots = knots, degree = 0,
                       Boundary.knots = c(0, 1)),
               "beyond boundary knots")
expect_warning(mSpline(c(x, 10), knots = knots, degree = 3,
                       Boundary.knots = c(0, 1)),
               "beyond boundary knots")
