## get implementations of v0.2.8 for reference
v2 <- new.env()
source("../v0.2.8.R", v2)
source("utils.R")

## helper functions
isNumMatrix <- v2$isNumMatrix

### 1. check correctness first
x <- seq.int(0, 1, 0.1)
knots <- c(0.3, 0.5, 0.6)
b_knots <- c(0, 1)

## without internal knots
nsMat0a <- nsk(x, df = 2, intercept = TRUE)
nsMat0b <- nsk(x, intercept = TRUE)

## integrals
nsMat1 <- nsk(x, intercept = FALSE, integral = TRUE)
## first derivatives
nsMat2 <- nsk(x, intercept = FALSE, derivs = 1)
## second derivatives
nsMat3 <- nsk(x, intercept = FALSE, derivs = 2)

## check matrix size
expect_true(isNumMatrix(nsMat0a, length(x), 2L))
expect_true(isNumMatrix(nsMat0b, length(x), 2L))
expect_true(isNumMatrix(nsMat1, length(x), 1L))
expect_true(isNumMatrix(nsMat2, length(x), 1L))
expect_true(isNumMatrix(nsMat3, length(x), 1L))

## natural spline basis
nsMat0 <- nsk(x, knots = knots, intercept = TRUE)
## integrals
nsMat1 <- nsk(x, knots = knots, intercept = TRUE, integral = TRUE)
## first derivatives
nsMat2 <- nsk(x, knots = knots, intercept = TRUE, derivs = 1)
## second derivatives
nsMat3 <- nsk(x, knots = knots, intercept = TRUE, derivs = 2)

## check matrix size
expect_true(isNumMatrix(nsMat0, length(x), length(knots) + 2L))
expect_true(isNumMatrix(nsMat1, length(x), length(knots) + 2L))
expect_true(isNumMatrix(nsMat2, length(x), length(knots) + 2L))
expect_true(isNumMatrix(nsMat3, length(x), length(knots) + 2L))

## update trim to 0
nsMat0 <- update(nsMat0, trim = 0)
nsMat1 <- update(nsMat1, trim = 0)
nsMat2 <- update(nsMat2, trim = 0)
nsMat3 <- update(nsMat3, trim = 0)

## check the returned H matrix in attribute
h_mat <- attr(nsMat0, "H")
bsMat <- bSpline(x, knots = knots, intercept = TRUE)
expect_eqt(nsMat0, bsMat %*% h_mat)
expect_eqt(nsMat1, ibs(x, knots = knots, intercept = TRUE) %*% h_mat)
expect_eqt(nsMat2, deriv(bsMat) %*% h_mat)
expect_eqt(nsMat3, deriv(bsMat, 2) %*% h_mat)

## specify df directly instead of knots
for (j in seq.int(2, 10)) {
    expect_true(isNumMatrix(
        nsk(x, df = j), length(x), j
    ))
}

## use the deriv method
expect_equivalent(nsMat0, deriv(nsMat1))
expect_equivalent(nsMat2, deriv(nsMat0))
expect_equivalent(nsMat3, deriv(nsMat2))
expect_equivalent(nsMat3, deriv(nsMat0, 2))

## check second derivatives at boundary knots
expect_true(all(abs(predict(nsMat3, 0)) < 1e-12))
expect_true(all(abs(predict(nsMat3, 1)) < 1e-12))

## keep names of x
names(x) <- sample(LETTERS, length(x), replace = TRUE)
expect_equal(rownames(nsk(x, df = 3)), names(x))

## for x outside of boundary
xx <- seq.int(- 1, 2, 0.05)
knots <- c(0.3, 0.4, 0.6, 0.8)
b_knots <- c(0, 1)

nsMat <- nsk(xx, knots = knots, intercept = TRUE,
             Boundary.knots = b_knots)
d1Mat <- nsk(xx, knots = knots, intercept = TRUE,
             Boundary.knots = b_knots, derivs = 1)
d2Mat <- nsk(xx, knots = knots, intercept = TRUE,
             Boundary.knots = b_knots, derivs = 2)
iMat <- nsk(xx, knots = knots, intercept = TRUE,
            Boundary.knots = b_knots, integral = TRUE)

expect_equal(nsMat[1:5, ] - nsMat[2:6, ], nsMat[2:6, ] - nsMat[3:7, ])
expect_equal(nsMat[seq.int(length(xx) - 10, length(xx) - 5), ] -
             nsMat[seq.int(length(xx) - 9, length(xx) - 4), ],
             nsMat[seq.int(length(xx) - 9, length(xx) - 4), ] -
             nsMat[seq.int(length(xx) - 8, length(xx) - 3), ])

expect_equal(d1Mat[1:5, ], d1Mat[2:6, ])
expect_equal(d1Mat[seq.int(length(xx) - 10, length(xx) - 5), ],
             d1Mat[seq.int(length(xx) - 9, length(xx) - 4), ])

expect_equivalent(d2Mat[1:5, ], matrix(0, ncol = ncol(d2Mat), nrow = 5))
expect_equivalent(tail(d2Mat, 5), matrix(0, ncol = ncol(d2Mat), nrow = 5))

expect_equivalent(iMat[1:5, ], matrix(0, ncol = ncol(iMat), nrow = 5))
## expect_equal(iMat[seq.int(length(xx) - 10, length(xx) - 5), 2] -
##              iMat[seq.int(length(xx) - 9, length(xx) - 4), 2],
##              iMat[seq.int(length(xx) - 9, length(xx) - 4), 2] -
##              iMat[seq.int(length(xx) - 8, length(xx) - 3), 2])
expect_true(all(
    iMat[seq.int(length(xx) - 10, length(xx) - 5), 5] -
    iMat[seq.int(length(xx) - 9, length(xx) - 4), 5] <
    iMat[seq.int(length(xx) - 9, length(xx) - 4), 5] -
    iMat[seq.int(length(xx) - 8, length(xx) - 3), 5]
))

### 2. checking inputs
x <- c(NA, seq.int(0, 0.5, 0.1), NA, seq.int(0.6, 1, 0.1), NA)
knots <- c(0.25, 0.5, 0.8)
x2 <- c(- 1, 2, x)
b_knots <- c(0, 1)

## expect errors
expect_error(nsk(x2, df = 5, derivs = - 1))
expect_error(nsk(x2, df = 1))
expect_error(nsk(rep(NA, 10), df = 2))

## make sure internal knots are unique
x1 <- c(rep(0, 100), runif(10))
expect_warning(nsk(x1, df = 5), "duplicated knots")
