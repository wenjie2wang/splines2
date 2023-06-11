## get implementations of v0.2.8 for reference
v2 <- new.env()
source("../v0.2.8.R", v2)
source("utils.R")

## helper functions
isNumMatrix <- v2$isNumMatrix

### 1. check correctness first
x <- c(NA, seq.int(0, 0.5, 0.1), NA, seq.int(0.6, 1, 0.1), NA)
knots <- c(0.25, 0.5, 0.75)
x2 <- c(- 1, 2, x)
b_knots <- c(0, 1)

## default cubic splines without internal knots
expect_eqt(bSpline(x), v2$bSpline(x))
expect_eqt(bSpline(x, derivs = 1), dbs(x))
expect_eqt(bSpline(x, derivs = 2), dbs(x, derivs = 2))
expect_eqt(bSpline(x, integral = TRUE), ibs(x))
expect_eqt(bSpline(x), bSpline(x, derivs = 1, integral = TRUE))

## cubic splines with specified df
expect_eqt(bSpline(x, df = 5),
           v2$bSpline(x, df = 5))

## cubic splines with specified internal knots
expect_eqt(bSpline(x, knots = knots),
           v2$bSpline(x, knots = knots))

## qudractic splines without internal knots
expect_eqt(bSpline(x, degree = 2L),
           v2$bSpline(x, degree = 2L))

## complete basis with intercept
expect_eqt(bSpline(x, intercept = TRUE),
           v2$bSpline(x, intercept = TRUE))

## specified knots
expect_eqt(bSpline(x, knots = knots, intercept = TRUE),
           v2$bSpline(x, knots = knots, intercept = TRUE))

## specified df
expect_eqt(bSpline(x, df = 6, intercept = TRUE),
           v2$bSpline(x, df = 6, intercept = TRUE))

## degree zero
knots2 <- seq.int(0.2, 0.8, 0.2)
expect_eqt(bSpline(x, knots = knots2, degree = 0),
           v2$bSpline(x, knots = knots2, degree = 0))
expect_eqt(bSpline(x, knots = knots2, degree = 0, intercept = TRUE),
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

## true closed-form formula given the all knots and degree
## test with two internal knots
x3 <- seq.int(0, 5, 0.1)
b0_1 <- function(x) as.numeric(x >= 0 & x < 1)
b0_2 <- function(x) as.numeric(x >= 1 & x < 3)
b0_3 <- function(x) as.numeric(x >= 3 & x <= 5)
expect_eqt(bSpline(x3, knots = c(1, 3), degree = 0L, intercept = TRUE),
           cbind(b0_1(x3), b0_2(x3), b0_3(x3)))

## x outside of boundary
suppressWarnings({
    expect_eqt(
        bSpline(x2, df = 6, degree = 3, Boundary.knots = b_knots),
        v2$bSpline(x2, df = 6, degree = 3, Boundary.knots = b_knots)
    )
})
suppressWarnings({
    expect_eqt(
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

options("splines2.warn.outside" = FALSE)
expect_silent(bSpline(c(x, 10), knots = knots, degree = 0,
                      Boundary.knots = c(0, 1)))
expect_silent(bSpline(c(x, 10), knots = knots, degree = 3,
                      Boundary.knots = c(0, 1)))
options("splines2.warn.outside" = TRUE)

## warning if it cannot set internal knots based on quantiles
expect_warning(bSpline(rep(0.5, 10), df = 10, Boundary.knots = c(0, 1)),
               pattern = "duplicated")
expect_warning(bSpline(c(0, rep(1, 10)), df = 4, Boundary.knots = c(0, 1)),
               pattern = "boundary")

### periodic B-splines
## with specified df
x <- c(seq.int(0, 3, 0.01), NA, seq.int(3, 4, 0.1), NA)
b_knots <- c(0, 1)
is_in <- function(x, a = 0, b = 1) {
    x >= a & x <= b
}
## without specified boundary knots
tmp <- bSpline(x, df = 8, degree = 2, periodic = TRUE)
expect_equal(attr(tmp, "Boundary.knots"), range(x, na.rm = TRUE))
expect_equal(length(attr(tmp, "knots")), 8)

## intercept = TRUE
## basis
res0 <- bSpline(x, df = 6, degree = 3, intercept = TRUE,
                Boundary.knots = b_knots, periodic = TRUE)
tmp <- bSpline(x, df = 6, degree = 3, intercept = TRUE,
               Boundary.knots = b_knots, periodic = TRUE,
               derivs = 1, integral = TRUE)
expect_equal(res0[1, ], res0[nrow(res0) - 1, ])
expect_true(isNumMatrix(res0, length(x), 6))
expect_eqt(matrix(predict(res0, 0.25), nrow = 4,
                  ncol = ncol(res0), byrow = TRUE),
           predict(res0, 0.25 + 1:4))
expect_true(all(is_in(attr(res0, "knots"), 0, 1)))
expect_eqt(res0, tmp)

## first derivatives
res1 <- bSpline(x, df = 6, degree = 3, intercept = TRUE,
                Boundary.knots = b_knots, periodic = TRUE, derivs = 1)
expect_equal(res1[1, ], res1[nrow(res1) - 1, ])
expect_true(isNumMatrix(res1, length(x), 6))
expect_true(all(is_in(attr(res1, "knots"), 0, 1)))
expect_eqt(deriv(res0), res1)
expect_eqt(deriv(tmp), res1)
expect_eqt(matrix(predict(res1, 0.25), nrow = 4,
                  ncol = ncol(res1), byrow = TRUE),
           predict(res1, 0.25 + 1:4))

## second derivatives
res2 <- bSpline(x, df = 6, degree = 3, intercept = TRUE,
                Boundary.knots = b_knots, periodic = TRUE, derivs = 2)
expect_equal(res2[1, ], res2[nrow(res2) - 1, ])
expect_true(isNumMatrix(res2, length(x), 6))
expect_true(all(is_in(attr(res2, "knots"), 0, 1)))
expect_eqt(deriv(res1), res2)
expect_eqt(deriv(res0, 2), res2)
expect_eqt(deriv(tmp, 2), res2)
expect_eqt(matrix(predict(res2, 0.25), nrow = 4,
                  ncol = ncol(res2), byrow = TRUE),
           predict(res2, 0.25 + 1:4))

## third derivatives
res3 <- bSpline(x, df = 6, degree = 3, intercept = TRUE,
                Boundary.knots = b_knots, periodic = TRUE, derivs = 3)
expect_equal(res3[1, ], res3[nrow(res3) - 1, ])
expect_true(isNumMatrix(res3, length(x), 6))
expect_true(all(is_in(attr(res3, "knots"), 0, 1)))
expect_eqt(deriv(res2), res3)
expect_eqt(deriv(res1, 2), res3)
expect_eqt(deriv(res0, 3), res3)
expect_eqt(deriv(tmp, 3), res3)
expect_eqt(matrix(predict(res3, 0.25), nrow = 4,
                  ncol = ncol(res3), byrow = TRUE),
           predict(res3, 0.25 + 1:4))

## fourth derivatives
res4 <- bSpline(x, df = 6, degree = 3, intercept = TRUE,
                Boundary.knots = b_knots, periodic = TRUE, derivs = 4)
expect_equal(res4[1, ], res4[nrow(res4) - 1, ])
expect_true(isNumMatrix(res4, length(x), 6))
expect_eqt(res4[1, , drop = FALSE],
           matrix(0, ncol = ncol(res4), nrow = 1))

## integrals
res3 <- bSpline(x, df = 6, degree = 3, intercept = TRUE,
                Boundary.knots = b_knots,
                periodic = TRUE, integral = TRUE)
expect_true(isNumMatrix(res3, length(x), 6))
expect_true(all(is_in(attr(res3, "knots"), 0, 1)))
expect_eqt(deriv(res3), res0)
expect_eqt(deriv(res3), tmp)
expect_eqt(deriv(res3, 2), res1)
expect_eqt(deriv(res3, 3), res2)
expect_eqt(matrix(predict(res3, 0.25), byrow = TRUE,
                  nrow = 4, ncol = ncol(res3)) +
           predict(res3, 1:4),
           predict(res3, 0.25 + 1:4))

## intercept = FALSE
## basis
res0 <- bSpline(x, df = 6, degree = 3, intercept = FALSE,
                Boundary.knots = b_knots, periodic = TRUE)
tmp <- bSpline(x, df = 6, degree = 3, intercept = FALSE,
               Boundary.knots = b_knots, periodic = TRUE,
               derivs = 1, integral = TRUE)
expect_equal(res0[1, ], res0[nrow(res0) - 1, ])
expect_eqt(matrix(predict(res0, 0.25), nrow = 4,
                  ncol = ncol(res0), byrow = TRUE),
           predict(res0, 0.25 + 1:4))
expect_eqt(res0, tmp)
expect_true(isNumMatrix(res0, length(x), 6))
expect_true(all(is_in(attr(res0, "knots"), 0, 1)))
## first derivatives
res1 <- bSpline(x, df = 6, degree = 3, intercept = FALSE,
                Boundary.knots = b_knots, periodic = TRUE, derivs = 1)
expect_equal(res1[1, ], res1[nrow(res1) - 1, ])
expect_true(isNumMatrix(res1, length(x), 6))
expect_true(all(is_in(attr(res1, "knots"), 0, 1)))
expect_eqt(deriv(res0), res1)
expect_eqt(deriv(tmp), res1)
expect_eqt(matrix(predict(res1, 0.25), nrow = 4,
                  ncol = ncol(res1), byrow = TRUE),
           predict(res1, 0.25 + 1:4))
## second derivatives
res2 <- bSpline(x, df = 6, degree = 3, intercept = FALSE,
                Boundary.knots = b_knots, periodic = TRUE, derivs = 2)
expect_equal(res2[1, ], res2[nrow(res2) - 1, ])
expect_true(isNumMatrix(res2, length(x), 6))
expect_true(all(is_in(attr(res2, "knots"), 0, 1)))
expect_eqt(deriv(res1), res2)
expect_eqt(deriv(res0, 2), res2)
expect_eqt(deriv(tmp, 2), res2)
expect_eqt(matrix(predict(res2, 0.25), nrow = 4,
                  ncol = ncol(res2), byrow = TRUE),
           predict(res2, 0.25 + 1:4))

## third derivatives
res3 <- bSpline(x, df = 6, degree = 3, intercept = FALSE,
                Boundary.knots = b_knots, periodic = TRUE, derivs = 3)
expect_equal(res3[1, ], res3[nrow(res3) - 1, ])
expect_true(isNumMatrix(res3, length(x), 6))
expect_true(all(is_in(attr(res3, "knots"), 0, 1)))
expect_eqt(deriv(res2), res3)
expect_eqt(deriv(res1, 2), res3)
expect_eqt(deriv(res0, 3), res3)
expect_eqt(deriv(tmp, 3), res3)
expect_eqt(matrix(predict(res3, 0.25), nrow = 4,
                  ncol = ncol(res3), byrow = TRUE),
           predict(res3, 0.25 + 1:4))

## fourth derivatives
res4 <- bSpline(x, df = 6, degree = 3, intercept = FALSE,
                Boundary.knots = b_knots, periodic = TRUE, derivs = 4)
expect_equal(res4[1, ], res4[nrow(res4) - 1, ])
expect_true(isNumMatrix(res4, length(x), 6))
expect_eqt(res4[1, , drop = FALSE],
           matrix(0, ncol = ncol(res4), nrow = 1))

## integrals
res3 <- bSpline(x, df = 6, degree = 3, intercept = FALSE,
                Boundary.knots = b_knots,
                periodic = TRUE, integral = TRUE)
expect_true(isNumMatrix(res3, length(x), 6))
expect_true(all(is_in(attr(res3, "knots"), 0, 1)))
expect_eqt(deriv(res3), res0)
expect_eqt(deriv(res3), tmp)
expect_eqt(deriv(res3, 2), res1)
expect_eqt(deriv(res3, 3), res2)
expect_eqt(matrix(predict(res3, 0.25), byrow = TRUE,
                  nrow = 4, ncol = ncol(res3)) +
           predict(res3, 1:4),
           predict(res3, 0.25 + 1:4))
