## get implementations of v0.2.8 for reference
v2 <- new.env()
source("../v0.2.8.R", v2)

## helper functions
isNumMatrix <- v2$isNumMatrix

### 1. check correctness first
## compare with splines::splineDesign
x <- seq.int(0, 1, 0.05)
ord <- 4
aKnots <- c(rep(0, ord), rep(1, ord))
expect_equivalent(dbs(x, derivs = 1, intercept = TRUE),
                  splines::splineDesign(aKnots, x = x, derivs = 1))
expect_equivalent(dbs(x, derivs = 2, intercept = TRUE),
                  splines::splineDesign(aKnots, x = x, derivs = 2))
## except at right boundary knots
expect_equivalent(
    dbs(x, derivs = 3, intercept = TRUE)[- length(x), ],
    splines::splineDesign(aKnots, derivs = 3, x = x)[- length(x), ]
)

knots <- c(0.2, 0.4, 0.7)
aKnots <- c(rep(0, ord), na.omit(knots), rep(1, ord))
expect_equivalent(dbs(x, derivs = 1, knots = knots, intercept = TRUE),
                  splines::splineDesign(aKnots, x = x, derivs = 1))
expect_equivalent(dbs(x, derivs = 2, knots = knots, intercept = TRUE),
                  splines::splineDesign(aKnots, x = x, derivs = 2))
expect_equivalent(
    dbs(x, derivs = 3, knots = knots, intercept = TRUE)[- length(x), ],
    splines::splineDesign(aKnots, x = x, derivs = 3)[- length(x), ]
)

knots <- c(0.3, 0.6)
ord <- 5
aKnots <- c(rep(0, ord), na.omit(knots), rep(1, ord))

expect_equivalent(dbs(x, 1, knots = knots, degree = 4, intercept = TRUE),
                  splines::splineDesign(aKnots, x, ord, derivs = 1))
expect_equivalent(dbs(x, 2, knots = knots, degree = 4, intercept = TRUE),
                  splines::splineDesign(aKnots, x, ord, derivs = 2))
expect_equivalent(dbs(x, 3, knots = knots, degree = 4, intercept = TRUE),
                  splines::splineDesign(aKnots, x, ord, derivs = 3))
expect_equivalent(
    dbs(x, 4, knots = knots, degree = 4, intercept = TRUE)[- length(x), ],
    splines::splineDesign(aKnots, x, ord, derivs = 4)[- length(x), ]
)

expect_error(dbs(x, 1, df = 1, intercept = TRUE))
expect_error(dbs(x, 1, df = 2, intercept = TRUE))
expect_error(dbs(x, 1, df = 3, intercept = TRUE))
expect_error(dbs(x, 2, df = 3, intercept = TRUE))
expect_true(
    isNumMatrix(dbs(x, 1, df = 1, degree = 0, intercept = TRUE), 21L, 1L)
)
expect_true(isNumMatrix(dbs(x, 1, df = 4), 21L, 4L))
expect_true(isNumMatrix(dbs(x, 1, df = 4, intercept = TRUE), 21L, 4L))
expect_true(isNumMatrix(dbs(x, 1, df = 5), 21L, 5L))
expect_true(isNumMatrix(dbs(x, 1, df = 5, intercept = TRUE), 21L, 5L))
expect_true(isNumMatrix(dbs(x, 1, df = 5, degree = 0), 21L, 5L))
expect_true(
    isNumMatrix(dbs(x, 1, df = 5, degree = 0, intercept = TRUE), 21L, 5L)
)

## compare with v0.2.8
x <- c(NA, seq.int(0, 0.5, 0.1), NA, seq.int(0.6, 1, 0.1), NA)
knots <- c(0.25, 0.5, 0.75)
x2 <- c(- 1, 2, x)
b_knots <- c(0, 1)

## default cubic splines without internal knots
expect_equivalent(dbs(x), v2$dbs(x))

## cubic splines with specified df
expect_equivalent(dbs(x, df = 5),
                  v2$dbs(x, df = 5))

## cubic splines with specified internal knots
expect_equivalent(dbs(x, knots = knots),
                  v2$dbs(x, knots = knots))

## qudractic splines without internal knots
expect_equivalent(dbs(x, degree = 2L),
                  v2$dbs(x, degree = 2L))

## complete basis with intercept
expect_equivalent(dbs(x, intercept = TRUE),
                  v2$dbs(x, intercept = TRUE))

## specified knots
expect_equivalent(dbs(x, knots = knots, intercept = TRUE),
                  v2$dbs(x, knots = knots, intercept = TRUE))

## specified df
expect_equivalent(dbs(x, df = 6, intercept = TRUE),
                  v2$dbs(x, df = 6, intercept = TRUE))

## degree zero
expect_equivalent(dbs(x, df = 5, degree = 0),
                  v2$dbs(x, df = 5, degree = 0))
expect_equivalent(dbs(x, df = 5, degree = 0, intercept = TRUE),
                  v2$dbs(x, df = 5, degree = 0, intercept = TRUE))
bsMat0a <- dbs(x, degree = 0, intercept = TRUE)
bsMat0b <- dbs(x, df = 5, degree = 0)
bsMat0c <- dbs(x, df = 5, degree = 0, intercept = TRUE)
bsMat0d <- dbs(x, knots = knots, degree = 0)
bsMat0e <- dbs(x, knots = knots, degree = 0, intercept = TRUE)
expect_true(isNumMatrix(bsMat0a, 14L, 1L))
expect_equal(sum(is.na(bsMat0b)), 15L) # keep NA's as is
expect_true(isNumMatrix(bsMat0b, 14L, 5L))
expect_true(isNumMatrix(bsMat0c, 14L, 5L))
expect_true(isNumMatrix(bsMat0d, 14L, 3L))
expect_true(isNumMatrix(bsMat0e, 14L, 4L))
expect_true(isNumMatrix(
    dbs(x, df = 10, knots = knots, degree = 0L),
    14L, 3L))
expect_true(isNumMatrix(
    dbs(x, df = 10, knots = knots,
        degree = 0, intercept = TRUE),
    14L, 4L))

## x outside of boundary
suppressWarnings({
    expect_equivalent(
        dbs(x2, df = 6, degree = 3, Boundary.knots = b_knots),
        v2$dbs(x2, df = 6, degree = 3, Boundary.knots = b_knots)
    )
})
suppressWarnings({
    expect_equivalent(
        dbs(x2, knots = knots, degree = 3, Boundary.knots = b_knots),
        v2$dbs(x2, knots = knots, degree = 3, Boundary.knots = b_knots)
    )
})

## keep names of x
names(x) <- sample(LETTERS, length(x), replace = TRUE)
expect_equal(rownames(dbs(x)), names(x))

### 2. check designed features with expectation
## NA is only allowed in x

## error if all of x are NA's
expect_error(dbs(c(NA_real_, NA_real_), degree = 0))
expect_error(dbs(c(NA, NA), df = 5))

## error if degree has NA or negative
expect_error(dbs(x, degree = - 1))
expect_error(dbs(x, degree = NA))

## error if df has NA or negative
expect_error(dbs(x, df = - 1))
expect_error(dbs(x, df = NA))

## error if knots has NA
expect_error(dbs(x, knots = c(0.1, 0.5, NA)))
expect_error(dbs(x, Boundary.knots = c(0.1, 0.5, NA)))

## error if boundary knots are inappropriate
expect_error(dbs(x, Boundary.knots = 0.1))
expect_error(dbs(x, Boundary.knots = c(0.1, 0.1)))
expect_error(dbs(x, Boundary.knots = c(0.1, 0.5, 1)))

## error if empty matrix
expect_true(isNumMatrix(dbs(x, degree = 0, intercept = TRUE),
                        length(x), 1))
expect_error(dbs(x, degree = 0))

## error if any internal knot is not placed inside of boundary
expect_error(dbs(x, knots = c(0, 0.5), degree = 0))

## warning if any x outside of boundary
expect_warning(dbs(c(x, 10), knots = knots, degree = 0,
                   Boundary.knots = c(0, 1)))
expect_warning(dbs(c(x, 10), knots = knots, degree = 3,
                   Boundary.knots = c(0, 1)))
