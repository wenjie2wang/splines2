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
expect_equivalent(cSpline(x), v2$cSpline(x))

## cubic splines with specified df
expect_equivalent(cSpline(x, df = 5),
                  v2$cSpline(x, df = 5))

## cubic splines with specified internal knots
expect_equivalent(cSpline(x, knots = knots),
                  v2$cSpline(x, knots = knots))

## qudractic splines without internal knots
expect_equivalent(cSpline(x, degree = 2L),
                  v2$cSpline(x, degree = 2L))

## complete basis with intercept
expect_equivalent(cSpline(x, intercept = TRUE),
                  v2$cSpline(x, intercept = TRUE))

## specified knots
expect_equivalent(cSpline(x, knots = knots, intercept = TRUE),
                  v2$cSpline(x, knots = knots, intercept = TRUE))

## specified df
expect_equivalent(cSpline(x, df = 6, intercept = TRUE),
                  v2$cSpline(x, df = 6, intercept = TRUE))

## degree zero
expect_equivalent(cSpline(x, df = 5, degree = 0),
                  v2$cSpline(x, df = 5, degree = 0))
expect_equivalent(cSpline(x, df = 5, degree = 0, intercept = TRUE),
                  v2$cSpline(x, df = 5, degree = 0, intercept = TRUE))
bsMat0a <- cSpline(x, degree = 0, intercept = TRUE)
bsMat0b <- cSpline(x, df = 5, degree = 0, intercept = FALSE)
bsMat0c <- cSpline(x, df = 5, degree = 0, intercept = TRUE)
bsMat0d <- cSpline(x, knots = knots, degree = 0, intercept = FALSE)
bsMat0e <- cSpline(x, knots = knots, degree = 0, intercept = TRUE)
expect_true(isNumMatrix(bsMat0a, 14L, 1L))
expect_equal(sum(is.na(bsMat0b)), 15L) # keep NA's as is
expect_true(isNumMatrix(bsMat0b, 14L, 5L))
expect_true(isNumMatrix(bsMat0c, 14L, 5L))
expect_true(isNumMatrix(bsMat0d, 14L, 3L))
expect_true(isNumMatrix(bsMat0e, 14L, 4L))
expect_true(isNumMatrix(
    cSpline(x, df = 10, knots = knots, degree = 0L, intercept = FALSE),
    14L, 3L))
expect_true(isNumMatrix(
    cSpline(x, df = 10, knots = knots,
            degree = 0, intercept = TRUE),
    14L, 4L))

## x outside of boundary
suppressWarnings({
    expect_equivalent(
        cSpline(x2, df = 6, degree = 3, Boundary.knots = b_knots),
        v2$cSpline(x2, df = 6, degree = 3, Boundary.knots = b_knots)
    )
})
suppressWarnings({
    expect_equivalent(
        cSpline(x2, knots = knots, degree = 3, Boundary.knots = b_knots),
        v2$cSpline(x2, knots = knots, degree = 3, Boundary.knots = b_knots)
    )
})

## keep names of x
names(x) <- sample(LETTERS, length(x), replace = TRUE)
expect_equal(rownames(cSpline(x)), names(x))

## equivalency with I-splines
expect_equivalent(
    cSpline(x, df = 5, derivs = 1,
            scale = FALSE, intercept = TRUE),
    iSpline(x, df = 5, intercept = TRUE)
)
expect_equivalent(
    cSpline(x, knots = knots, derivs = 1,
            scale = FALSE, intercept = TRUE),
    iSpline(x, knots = knots, intercept = TRUE)
)
expect_equivalent(
    cSpline(x, df = 6, degree = 2, derivs = 2,
            intercept = TRUE, scale = FALSE),
    iSpline(x, df = 6, degree = 2, derivs = 1,
            intercept = TRUE)
)
expect_equivalent(
    cSpline(x, knots = knots, degree = 2, derivs = 2,
            intercept = TRUE, scale = FALSE),
    iSpline(x, knots = knots, degree = 2, derivs = 1,
            intercept = TRUE)
)
expect_equivalent(
    cSpline(x, df = 6, degree = 2, derivs = 3,
            intercept = TRUE, scale = FALSE),
    iSpline(x, df = 6, degree = 2, derivs = 2,
            intercept = TRUE)
)
expect_equivalent(
    cSpline(x, knots = knots, degree = 2, derivs = 3,
            intercept = TRUE, scale = FALSE),
    iSpline(x, knots = knots, degree = 2, derivs = 2,
            intercept = TRUE)
)
expect_equivalent(
    cSpline(x, df = 6, degree = 2, derivs = 4,
            intercept = TRUE, scale = FALSE),
    iSpline(x, df = 6, degree = 2, derivs = 3,
            intercept = TRUE)
)
expect_equivalent(
    cSpline(x, knots = knots, degree = 2, derivs = 4,
            intercept = TRUE, scale = FALSE),
    iSpline(x, knots = knots, degree = 2, derivs = 3,
            intercept = TRUE)
)


### 2. check designed features with expectation
## NA is only allowed in x

## error if all of x are NA's
expect_error(cSpline(c(NA_real_, NA_real_), degree = 0))
expect_error(cSpline(c(NA, NA), df = 5))

## error if degree has NA or negative
expect_error(cSpline(x, degree = - 1))
expect_error(cSpline(x, degree = NA))

## error if df has NA or negative
expect_error(cSpline(x, df = - 1))
expect_error(cSpline(x, df = NA))

## error if knots has NA
expect_error(cSpline(x, knots = c(0.1, 0.5, NA)))
expect_error(cSpline(x, Boundary.knots = c(0.1, 0.5, NA)))

## error if boundary knots are inappropriate
expect_error(cSpline(x, Boundary.knots = 0.1))
expect_error(cSpline(x, Boundary.knots = c(0.1, 0.1)))
expect_error(cSpline(x, Boundary.knots = c(0.1, 0.5, 1)))

## error if empty matrix
expect_true(isNumMatrix(cSpline(x, degree = 0, intercept = TRUE),
                        length(x), 1))
expect_error(cSpline(x, degree = 0, intercept = FALSE))

## error if any internal knot is not placed inside of boundary
expect_error(cSpline(x, knots = c(0, 0.5), degree = 0))

## warning if any x outside of boundary
expect_warning(cSpline(c(x, 10), knots = knots, degree = 0,
                       Boundary.knots = c(0, 1)))
expect_warning(cSpline(c(x, 10), knots = knots, degree = 3,
                       Boundary.knots = c(0, 1)))
