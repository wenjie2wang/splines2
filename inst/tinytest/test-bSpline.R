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
expect_equivalent(bSpline(x), v2$bSpline(x))

## cubic splines with specified df
expect_equivalent(bSpline(x, df = 5),
                  v2$bSpline(x, df = 5))

## cubic splines with specified internal knots
expect_equivalent(bSpline(x, knots = knots),
                  v2$bSpline(x, knots = knots))

## qudractic splines without internal knots
expect_equivalent(bSpline(x, degree = 2L),
                  v2$bSpline(x, degree = 2L))

## complete basis with intercept
expect_equivalent(bSpline(x, intercept = TRUE),
                  v2$bSpline(x, intercept = TRUE))

## specified knots
expect_equivalent(bSpline(x, knots = knots, intercept = TRUE),
                  v2$bSpline(x, knots = knots, intercept = TRUE))

## specified df
expect_equivalent(bSpline(x, df = 6, intercept = TRUE),
                  v2$bSpline(x, df = 6, intercept = TRUE))

## degree zero
knots2 <- seq.int(0.2, 0.8, 0.2)
expect_equivalent(bSpline(x, knots = knots2, degree = 0),
                  v2$bSpline(x, knots = knots2, degree = 0))
expect_equivalent(bSpline(x, knots = knots2, degree = 0, intercept = TRUE),
                  v2$bSpline(x, knots = knots2, degree = 0, intercept = TRUE))
bsMat0a <- bSpline(x, degree = 0, intercept = TRUE)
bsMat0b <- bSpline(x, df = 5, degree = 0)
bsMat0c <- bSpline(x, df = 5, degree = 0, intercept = TRUE)
bsMat0d <- bSpline(x, knots = knots, degree = 0)
bsMat0e <- bSpline(x, knots = knots, degree = 0, intercept = TRUE)
expect_true(isNumMatrix(bsMat0a, 14L, 1L))
expect_equal(sum(is.na(bsMat0b)), 15L) # keep NA's as is
expect_true(isNumMatrix(bsMat0b, 14L, 5L))
expect_true(isNumMatrix(bsMat0c, 14L, 5L))
expect_true(isNumMatrix(bsMat0d, 14L, 3L))
expect_true(isNumMatrix(bsMat0e, 14L, 4L))
expect_true(isNumMatrix(
    bSpline(x, df = 10, knots = knots, degree = 0L),
    14L, 3L))
expect_true(isNumMatrix(
    bSpline(x, df = 10, knots = knots,
            degree = 0, intercept = TRUE),
    14L, 4L))

## true close form formula given the all knots and degree
## test with two internal knots
x3 <- seq.int(0, 5, 0.1)
b0_1 <- function(x) as.numeric(x >= 0 & x < 1)
b0_2 <- function(x) as.numeric(x >= 1 & x < 3)
b0_3 <- function(x) as.numeric(x >= 3 & x <= 5)
expect_equivalent(bSpline(x3, knots = c(1, 3), degree = 0L, intercept = TRUE),
                  cbind(b0_1(x3), b0_2(x3), b0_3(x3)))

## x outside of boundary
suppressWarnings({
    expect_equivalent(
        bSpline(x2, df = 6, degree = 3, Boundary.knots = b_knots),
        v2$bSpline(x2, df = 6, degree = 3, Boundary.knots = b_knots)
    )
})
suppressWarnings({
    expect_equivalent(
        bSpline(x2, knots = knots, degree = 3, Boundary.knots = b_knots),
        v2$bSpline(x2, knots = knots, degree = 3, Boundary.knots = b_knots)
    )
})

## keep names of x
names(x) <- sample(LETTERS, length(x), replace = TRUE)
expect_equal(rownames(bSpline(x)), names(x))

### 2. check designed features with expectation
## NA is only allowed in x

## error if all of x are NA's
expect_error(bSpline(c(NA_real_, NA_real_), degree = 0))
expect_error(bSpline(c(NA, NA), df = 5))

## error if degree has NA or negative
expect_error(bSpline(x, degree = - 1))
expect_error(bSpline(x, degree = NA))

## error if df has NA or negative
expect_error(bSpline(x, df = - 1))
expect_error(bSpline(x, df = NA))

## error if knots has NA
expect_error(bSpline(x, knots = c(0.1, 0.5, NA)))
expect_error(bSpline(x, Boundary.knots = c(0.1, 0.5, NA)))

## error if boundary knots are inappropriate
expect_error(bSpline(x, Boundary.knots = 0.1))
expect_error(bSpline(x, Boundary.knots = c(0.1, 0.1)))
expect_error(bSpline(x, Boundary.knots = c(0.1, 0.5, 1)))

## error if empty matrix
expect_true(isNumMatrix(bSpline(x, degree = 0, intercept = TRUE),
                        length(x), 1))
expect_error(bSpline(x, degree = 0))

## error if any internal knot is placed outside boundary
expect_error(bSpline(x, knots = c(- 0.1, 0.5), degree = 0))
expect_error(bSpline(x, knots = c(range(x), 0.5), degree = 0))

## warning if any x outside of boundary
expect_warning(bSpline(c(x, 10), knots = knots, degree = 0,
                       Boundary.knots = c(0, 1)))
expect_warning(bSpline(c(x, 10), knots = knots, degree = 3,
                       Boundary.knots = c(0, 1)))
