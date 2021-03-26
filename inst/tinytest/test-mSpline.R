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
msMat0b2 <- mSpline(x, knots = knots, degree = 1, derivs = 1, integral = TRUE)
isMat0a <- mSpline(x, knots = knots, degree = 1, integral = TRUE)
msMat0c <- mSpline(x, knots = knots, degree = 1, intercept = TRUE)
msMat0d <- mSpline(x, knots = knots, degree = 2)
msMat0e <- mSpline(x, knots = knots, degree = 2, intercept = TRUE)
msMat0f <- mSpline(0.1, knots = knots, degree = 2, intercept = TRUE,
                   Boundary.knots = c(0, 1), derivs = 1)
msMat0g <- mSpline(x, df = 6, degree = 2, intercept = TRUE,
                   Boundary.knots = c(0, 1), derivs = 2)
msMat0h <- mSpline(0.1, knots = knots, degree = 2,
                   Boundary.knots = c(0, 1), derivs = 3)
expect_true(isNumMatrix(msMat0a, 14L, 1L))
expect_equal(sum(is.na(msMat0b)), 12L) # keep NA's as is
expect_equivalent(msMat0b, msMat0b2)
expect_equivalent(deriv(isMat0a), msMat0b2)
expect_true(isNumMatrix(msMat0b, 14L, 4L))
expect_true(isNumMatrix(msMat0c, 14L, 5L))
expect_true(isNumMatrix(msMat0d, 14L, 5L))
expect_true(isNumMatrix(msMat0e, 14L, 6L))
expect_true(isNumMatrix(msMat0f, 1L, 6L))
expect_true(isNumMatrix(msMat0g, 14L, 6L))
expect_true(isNumMatrix(msMat0h, 1L, 5L))
expect_error(mSpline(x, degree = 0))
expect_warning(mSpline(c(x, 10), knots = knots, degree = 0,
                       Boundary.knots = c(0, 1)),
               "beyond boundary knots")

## intercept = FALSE
msMat0f <- mSpline(0.1, knots = knots, degree = 2,
                   Boundary.knots = c(0, 1), derivs = 1)
msMat0g <- mSpline(x, df = 6, degree = 2,
                   Boundary.knots = c(0, 1), derivs = 2)
expect_true(isNumMatrix(msMat0f, 1, length(knots) + 2))
expect_true(isNumMatrix(msMat0g, length(x), 6))

## true close form formula given the all knots and degree
## transformation of constant basis
x3 <- seq.int(0, 7, 0.1)
m0_1 <- function(x) as.numeric(x < 1)
m0_2 <- function(x) as.numeric(x >= 1 & x < 3) * 0.5
m0_3 <- function(x) as.numeric(x >= 3 & x <= 7) * 0.25
expect_equivalent(mSpline(x3, knots = c(1, 3), degree = 0L, intercept = TRUE),
                  cbind(m0_1(x3), m0_2(x3), m0_3(x3)))

## transformation of linear basis
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
expect_error(mSpline(x, knots = c(0.1, 0.5, NA)))
expect_error(mSpline(x, Boundary.knots = c(0.1, 0.5, NA)))

## error if boundary knots are inappropriate
expect_error(mSpline(x, Boundary.knots = 0.1))
expect_error(mSpline(x, Boundary.knots = c(0.1, 0.1)))
expect_error(mSpline(x, Boundary.knots = c(0.1, 0.5, 1)))

## error if empty matrix
expect_true(isNumMatrix(mSpline(x, degree = 0, intercept = TRUE),
                        length(x), 1))
expect_error(mSpline(x, degree = 0))
expect_error(mSpline(x, degree = 0, derivs = 1))

## error if any internal knot is placed outside boundary
expect_error(mSpline(x, knots = c(- 0.1, 0.5), degree = 0))

## warning if any x outside of boundary
expect_warning(mSpline(c(x, 10), knots = knots, degree = 0,
                       Boundary.knots = c(0, 1)))
expect_warning(mSpline(c(x, 10), knots = knots, degree = 3,
                       Boundary.knots = c(0, 1)))

### 3. periodic M-splines
## 3.1. with specified df
x <- c(seq.int(0, 3, 0.01), NA, seq.int(3, 4, 0.1), NA)
b_knots <- c(0, 1)

is_in <- function(x, a = 0, b = 1) {
    x >= a & x <= b
}

## without specified boundary knots
tmp <- mSpline(x, df = 8, degree = 2, periodic = TRUE)
expect_equal(attr(tmp, "Boundary.knots"), range(x, na.rm = TRUE))
expect_equal(length(attr(tmp, "knots")), 8)

## intercept = TRUE
## basis
res0 <- mSpline(x, df = 6, degree = 3, intercept = TRUE,
                Boundary.knots = b_knots, periodic = TRUE)
tmp <- mSpline(x, df = 6, degree = 3, intercept = TRUE,
               Boundary.knots = b_knots, periodic = TRUE,
               derivs = 1, integral = TRUE)
expect_equal(res0[1, ], res0[nrow(res0) - 1, ])
expect_true(isNumMatrix(res0, length(x), 6))
expect_equivalent(matrix(predict(res0, 0.25), nrow = 4,
                         ncol = ncol(res0), byrow = TRUE),
                  predict(res0, 0.25 + 1:4))
expect_true(all(is_in(attr(res0, "knots"), 0, 1)))
expect_equivalent(res0, tmp)

## first derivatives
res1 <- mSpline(x, df = 6, degree = 3, intercept = TRUE,
                Boundary.knots = b_knots, periodic = TRUE, derivs = 1)
expect_equal(res1[1, ], res1[nrow(res1) - 1, ])
expect_true(isNumMatrix(res1, length(x), 6))
expect_true(all(is_in(attr(res1, "knots"), 0, 1)))
expect_equivalent(deriv(res0), res1)
expect_equivalent(deriv(tmp), res1)
expect_equivalent(matrix(predict(res1, 0.25), nrow = 4,
                         ncol = ncol(res1), byrow = TRUE),
                  predict(res1, 0.25 + 1:4))

## second derivatives
res2 <- mSpline(x, df = 6, degree = 3, intercept = TRUE,
                Boundary.knots = b_knots, periodic = TRUE, derivs = 2)
expect_equal(res2[1, ], res2[nrow(res2) - 1, ])
expect_true(isNumMatrix(res2, length(x), 6))
expect_true(all(is_in(attr(res2, "knots"), 0, 1)))
expect_equivalent(deriv(res1), res2)
expect_equivalent(deriv(res0, 2), res2)
expect_equivalent(deriv(tmp, 2), res2)
expect_equivalent(matrix(predict(res2, 0.25), nrow = 4,
                         ncol = ncol(res2), byrow = TRUE),
                  predict(res2, 0.25 + 1:4))

## third derivatives
res3 <- mSpline(x, df = 6, degree = 3, intercept = TRUE,
                Boundary.knots = b_knots, periodic = TRUE, derivs = 3)
expect_equal(res3[1, ], res3[nrow(res3) - 1, ])
expect_true(isNumMatrix(res3, length(x), 6))
expect_true(all(is_in(attr(res3, "knots"), 0, 1)))
expect_equivalent(deriv(res2), res3)
expect_equivalent(deriv(res1, 2), res3)
expect_equivalent(deriv(res0, 3), res3)
expect_equivalent(deriv(tmp, 3), res3)
expect_equivalent(matrix(predict(res3, 0.25), nrow = 4,
                         ncol = ncol(res3), byrow = TRUE),
                  predict(res3, 0.25 + 1:4))

## fourth derivatives
res4 <- mSpline(x, df = 6, degree = 3, intercept = TRUE,
                Boundary.knots = b_knots, periodic = TRUE, derivs = 4)
expect_equal(res4[1, ], res4[nrow(res4) - 1, ])
expect_true(isNumMatrix(res4, length(x), 6))
expect_equivalent(res4[1, , drop = FALSE],
                  matrix(0, ncol = ncol(res4), nrow = 1))

## integrals
res3 <- mSpline(x, df = 6, degree = 3, intercept = TRUE,
                Boundary.knots = b_knots,
                periodic = TRUE, integral = TRUE)
expect_true(isNumMatrix(res3, length(x), 6))
expect_true(all(is_in(attr(res3, "knots"), 0, 1)))
expect_equivalent(deriv(res3), res0)
expect_equivalent(deriv(res3), tmp)
expect_equivalent(deriv(res3, 2), res1)
expect_equivalent(deriv(res3, 3), res2)
expect_equivalent(matrix(predict(res3, 0.25), byrow = TRUE,
                         nrow = 4, ncol = ncol(res3)) + seq_len(4),
                  predict(res3, 0.25 + 1:4))

## intercept = FALSE
## basis
res0 <- mSpline(x, df = 6, degree = 3, intercept = FALSE,
                Boundary.knots = b_knots, periodic = TRUE)
tmp <- mSpline(x, df = 6, degree = 3, intercept = FALSE,
               Boundary.knots = b_knots, periodic = TRUE,
               derivs = 1, integral = TRUE)
expect_equal(res0[1, ], res0[nrow(res0) - 1, ])
expect_equivalent(matrix(predict(res0, 0.25), nrow = 4,
                         ncol = ncol(res0), byrow = TRUE),
                  predict(res0, 0.25 + 1:4))
expect_equivalent(res0, tmp)
expect_true(isNumMatrix(res0, length(x), 6))
expect_true(all(is_in(attr(res0, "knots"), 0, 1)))
## first derivatives
res1 <- mSpline(x, df = 6, degree = 3, intercept = FALSE,
                Boundary.knots = b_knots, periodic = TRUE, derivs = 1)
expect_equal(res1[1, ], res1[nrow(res1) - 1, ])
expect_true(isNumMatrix(res1, length(x), 6))
expect_true(all(is_in(attr(res1, "knots"), 0, 1)))
expect_equivalent(deriv(res0), res1)
expect_equivalent(deriv(tmp), res1)
expect_equivalent(matrix(predict(res1, 0.25), nrow = 4,
                         ncol = ncol(res1), byrow = TRUE),
                  predict(res1, 0.25 + 1:4))
## second derivatives
res2 <- mSpline(x, df = 6, degree = 3, intercept = FALSE,
                Boundary.knots = b_knots, periodic = TRUE, derivs = 2)
expect_equal(res2[1, ], res2[nrow(res2) - 1, ])
expect_true(isNumMatrix(res2, length(x), 6))
expect_true(all(is_in(attr(res2, "knots"), 0, 1)))
expect_equivalent(deriv(res1), res2)
expect_equivalent(deriv(res0, 2), res2)
expect_equivalent(deriv(tmp, 2), res2)
expect_equivalent(matrix(predict(res2, 0.25), nrow = 4,
                         ncol = ncol(res2), byrow = TRUE),
                  predict(res2, 0.25 + 1:4))

## third derivatives
res3 <- mSpline(x, df = 6, degree = 3, intercept = FALSE,
                Boundary.knots = b_knots, periodic = TRUE, derivs = 3)
expect_equal(res3[1, ], res3[nrow(res3) - 1, ])
expect_true(isNumMatrix(res3, length(x), 6))
expect_true(all(is_in(attr(res3, "knots"), 0, 1)))
expect_equivalent(deriv(res2), res3)
expect_equivalent(deriv(res1, 2), res3)
expect_equivalent(deriv(res0, 3), res3)
expect_equivalent(deriv(tmp, 3), res3)
expect_equivalent(matrix(predict(res3, 0.25), nrow = 4,
                         ncol = ncol(res3), byrow = TRUE),
                  predict(res3, 0.25 + 1:4))

## fourth derivatives
res4 <- mSpline(x, df = 6, degree = 3, intercept = FALSE,
                Boundary.knots = b_knots, periodic = TRUE, derivs = 4)
expect_equal(res4[1, ], res4[nrow(res4) - 1, ])
expect_true(isNumMatrix(res4, length(x), 6))
expect_equivalent(res4[1, , drop = FALSE],
                  matrix(0, ncol = ncol(res4), nrow = 1))

## integrals
res3 <- mSpline(x, df = 6, degree = 3, intercept = FALSE,
                Boundary.knots = b_knots,
                periodic = TRUE, integral = TRUE)
expect_true(isNumMatrix(res3, length(x), 6))
expect_true(all(is_in(attr(res3, "knots"), 0, 1)))
expect_equivalent(deriv(res3), res0)
expect_equivalent(deriv(res3), tmp)
expect_equivalent(deriv(res3, 2), res1)
expect_equivalent(deriv(res3, 3), res2)
expect_equivalent(matrix(predict(res3, 0.25), byrow = TRUE,
                         nrow = 4, ncol = ncol(res3)) + seq_len(4),
                  predict(res3, 0.25 + 1:4))

## 3.2. with specified knots
x <- c(seq.int(0, 3, 0.01), NA, seq.int(3, 4, 0.1), NA)
knots <- c(0.3, 0.6, 0.8)
b_knots <- c(0, 1)

## intercept = TRUE
## basis
res0 <- mSpline(x, knots = knots, degree = 3, intercept = TRUE,
                Boundary.knots = b_knots, periodic = TRUE)
tmp <- mSpline(x, knots = knots, degree = 3, intercept = TRUE,
               Boundary.knots = b_knots, periodic = TRUE,
               derivs = 1, integral = TRUE)
expect_equal(res0[1, ], res0[nrow(res0) - 1, ])
expect_equivalent(res0, tmp)
expect_true(isNumMatrix(res0, length(x), length(knots) + 1L))
expect_equivalent(matrix(predict(res0, 0.25), nrow = 4,
                         ncol = ncol(res0), byrow = TRUE),
                  predict(res0, 0.25 + 1:4))

## first derivatives
res1 <- mSpline(x, knots = knots, degree = 3, intercept = TRUE,
                Boundary.knots = b_knots, periodic = TRUE, derivs = 1)
expect_equal(res1[1, ], res1[nrow(res1) - 1, ])
expect_true(isNumMatrix(res1, length(x), length(knots) + 1L))
expect_equivalent(deriv(res0), res1)
expect_equivalent(deriv(tmp), res1)
expect_equivalent(matrix(predict(res1, 0.25), nrow = 4,
                         ncol = ncol(res1), byrow = TRUE),
                  predict(res1, 0.25 + 1:4))

## second derivatives
res2 <- mSpline(x, knots = knots, degree = 3, intercept = TRUE,
                Boundary.knots = b_knots, periodic = TRUE, derivs = 2)
expect_equal(res2[1, ], res2[nrow(res2) - 1, ])
expect_true(isNumMatrix(res2, length(x), length(knots) + 1L))
expect_equivalent(deriv(res1), res2)
expect_equivalent(deriv(res0, 2), res2)
expect_equivalent(deriv(tmp, 2), res2)
expect_equivalent(matrix(predict(res2, 0.25), nrow = 4,
                         ncol = ncol(res2), byrow = TRUE),
                  predict(res2, 0.25 + 1:4))

## integrals
res3 <- mSpline(x, knots = knots, degree = 3, intercept = TRUE,
                Boundary.knots = b_knots,
                periodic = TRUE, integral = TRUE)
expect_true(isNumMatrix(res3, length(x), length(knots) + 1L))
expect_equivalent(deriv(res3), res0)
expect_equivalent(deriv(res3), tmp)
expect_equivalent(deriv(res3, 2), res1)
expect_equivalent(deriv(res3, 3), res2)
expect_equivalent(matrix(predict(res3, 0.25), byrow = TRUE,
                         nrow = 4, ncol = ncol(res3)) + seq_len(4),
                  predict(res3, 0.25 + 1:4))

## intercept = FALSE
knots <- c(0.2, 0.5, 0.6, 0.75)
## basis
res0 <- mSpline(x, knots = knots, degree = 3, intercept = FALSE,
                Boundary.knots = b_knots, periodic = TRUE)
tmp <- mSpline(x, knots = knots, degree = 3, intercept = FALSE,
               Boundary.knots = b_knots, periodic = TRUE,
               derivs = 1, integral = TRUE)
expect_true(isNumMatrix(res0, length(x), length(knots)))
expect_equivalent(res0, tmp)
expect_equivalent(matrix(predict(res0, 0.25), nrow = 4,
                         ncol = ncol(res0), byrow = TRUE),
                  predict(res0, 0.25 + 1:4))

## first derivatives
res1 <- mSpline(x, knots = knots, degree = 3, intercept = FALSE,
                Boundary.knots = b_knots, periodic = TRUE, derivs = 1)
expect_equal(res1[1, ], res1[nrow(res1) - 1, ])
expect_true(isNumMatrix(res1, length(x), length(knots)))
expect_equivalent(deriv(res0), res1)
expect_equivalent(deriv(tmp), res1)
expect_equivalent(matrix(predict(res1, 0.25), nrow = 4,
                         ncol = ncol(res1), byrow = TRUE),
                  predict(res1, 0.25 + 1:4))

## second derivatives
res2 <- mSpline(x, knots = knots, degree = 3, intercept = FALSE,
                Boundary.knots = b_knots, periodic = TRUE, derivs = 2)
expect_equal(res2[1, ], res2[nrow(res2) - 1, ])
expect_true(isNumMatrix(res2, length(x), length(knots)))
expect_equivalent(deriv(res1), res2)
expect_equivalent(deriv(res0, 2), res2)
expect_equivalent(deriv(tmp, 2), res2)
expect_equivalent(matrix(predict(res2, 0.25), nrow = 4,
                         ncol = ncol(res2), byrow = TRUE),
                  predict(res2, 0.25 + 1:4))

## integrals
res3 <- mSpline(x, knots = knots, degree = 3, intercept = FALSE,
                Boundary.knots = b_knots,
                periodic = TRUE, integral = TRUE)
expect_true(isNumMatrix(res3, length(x), length(knots)))
expect_equivalent(deriv(res3), res0)
expect_equivalent(deriv(res3), tmp)
expect_equivalent(deriv(res3, 2), res1)
expect_equivalent(deriv(res3, 3), res2)
expect_equivalent(matrix(predict(res3, 0.25), byrow = TRUE,
                         nrow = 4, ncol = ncol(res3)) + seq_len(4),
                  predict(res3, 0.25 + 1:4))

### 4. catch errors and warnings for periodic splines
expect_error(mSpline(x, df = 2, periodic = TRUE))
expect_error(mSpline(x, df = 1, degree = 2, periodic = TRUE))
