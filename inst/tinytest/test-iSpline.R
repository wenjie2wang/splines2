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
expect_equivalent(iSpline(x), v2$iSpline(x))

## cubic splines with specified df
expect_equivalent(iSpline(x, df = 5),
                  v2$iSpline(x, df = 5))

## cubic splines with specified internal knots
expect_equivalent(iSpline(x, knots = knots),
                  v2$iSpline(x, knots = knots))

## qudractic splines without internal knots
expect_equivalent(iSpline(x, degree = 2L),
                  v2$iSpline(x, degree = 2L))

## complete basis with intercept
expect_equivalent(iSpline(x, intercept = TRUE),
                  v2$iSpline(x, intercept = TRUE))

## specified knots
expect_equivalent(iSpline(x, knots = knots, intercept = TRUE),
                  v2$iSpline(x, knots = knots, intercept = TRUE))

## specified df
expect_equivalent(iSpline(x, df = 6, intercept = TRUE),
                  v2$iSpline(x, df = 6, intercept = TRUE))

## degree zero
expect_equivalent(iSpline(x, df = 5, degree = 0),
                  v2$iSpline(x, df = 5, degree = 0))
expect_equivalent(iSpline(x, df = 5, degree = 0, intercept = TRUE),
                  v2$iSpline(x, df = 5, degree = 0, intercept = TRUE))
bsMat0a <- iSpline(x, degree = 0, intercept = TRUE)
bsMat0b <- iSpline(x, df = 5, degree = 0, intercept = FALSE)
bsMat0c <- iSpline(x, df = 5, degree = 0, intercept = TRUE)
bsMat0d <- iSpline(x, knots = knots, degree = 0, intercept = FALSE)
bsMat0e <- iSpline(x, knots = knots, degree = 0, intercept = TRUE)
expect_true(isNumMatrix(bsMat0a, 14L, 1L))
expect_equal(sum(is.na(bsMat0b)), 15L) # keep NA's as is
expect_true(isNumMatrix(bsMat0b, 14L, 5L))
expect_true(isNumMatrix(bsMat0c, 14L, 5L))
expect_true(isNumMatrix(bsMat0d, 14L, 3L))
expect_true(isNumMatrix(bsMat0e, 14L, 4L))
expect_true(isNumMatrix(
    iSpline(x, df = 10, knots = knots, degree = 0L, intercept = FALSE),
    14L, 3L))
expect_true(isNumMatrix(
    iSpline(x, df = 10, knots = knots,
            degree = 0, intercept = TRUE),
    14L, 4L))

## x outside of boundary
suppressWarnings({
    expect_equivalent(
        iSpline(x2, df = 6, degree = 3, Boundary.knots = b_knots),
        v2$iSpline(x2, df = 6, degree = 3, Boundary.knots = b_knots)
    )
})
suppressWarnings({
    expect_equivalent(
        iSpline(x2, knots = knots, degree = 3, Boundary.knots = b_knots),
        v2$iSpline(x2, knots = knots, degree = 3, Boundary.knots = b_knots)
    )
})

## keep names of x
names(x) <- sample(LETTERS, length(x), replace = TRUE)
expect_equal(rownames(iSpline(x)), names(x))

## equivalency with M-splines
expect_equivalent(
    iSpline(x, df = 5, derivs = 1),
    mSpline(x, df = 5, intercept = TRUE)
)
expect_equivalent(
    iSpline(x, knots = knots, degree = 2, derivs = 2),
    mSpline(x, knots = knots, degree = 2, derivs = 1, intercept = TRUE)
)
expect_equivalent(
    iSpline(x, knots = knots, degree = 2, derivs = 3),
    mSpline(x, knots = knots, degree = 2, derivs = 2, intercept = TRUE)
)


### 2. check designed features with expectation
## NA is only allowed in x

## error if all of x are NA's
expect_error(iSpline(c(NA_real_, NA_real_), degree = 0))
expect_error(iSpline(c(NA, NA), df = 5))

## error if degree has NA or negative
expect_error(iSpline(x, degree = - 1))
expect_error(iSpline(x, degree = NA))

## error if df has NA or negative
expect_error(iSpline(x, df = - 1))
expect_error(iSpline(x, df = NA))

## error if knots has NA
expect_error(iSpline(x, knots = c(0.1, 0.5, NA)))
expect_error(iSpline(x, Boundary.knots = c(0.1, 0.5, NA)))

## error if boundary knots are inappropriate
expect_error(iSpline(x, Boundary.knots = 0.1))
expect_error(iSpline(x, Boundary.knots = c(0.1, 0.1)))
expect_error(iSpline(x, Boundary.knots = c(0.1, 0.5, 1)))

## error if empty matrix
expect_true(isNumMatrix(iSpline(x, degree = 0, intercept = TRUE),
                        length(x), 1))
expect_error(iSpline(x, degree = 0, intercept = FALSE))

## error if any internal knot is placed outside boundary
expect_error(iSpline(x, knots = c(- 0.1, 0.5), degree = 0))

## warning if any x outside of boundary
expect_warning(iSpline(c(x, 10), knots = knots, degree = 0,
                       Boundary.knots = c(0, 1)))
expect_warning(iSpline(c(x, 10), knots = knots, degree = 3,
                       Boundary.knots = c(0, 1)))

## error for negative derivs
expect_error(iSpline(x, derivs = - 1))
