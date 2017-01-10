context("Testing bSpline")


test_that("call splines::bs", {
    x <- c(NA, seq.int(0, 0.5, 0.1), NA, seq.int(0.6, 1, 0.1), NA)
    knots <- c(0.25, 0.5, 0.75)
    ## for ease of testing
    bsFun <- function(x, df, knots, degree, intercept, Boundary.knots) {
        funCall <- match.call()
        funCall[[1L]] <- quote(splines::bs)
        bsMat <- eval(funCall)
        class(bsMat) <- "matrix"
        bsMat
    }
    expect_equivalent(bSpline(x), bsFun(x))
    expect_equivalent(bSpline(x, df = 5), bsFun(x, df = 5))
    expect_equivalent(bSpline(x, knots = knots), bsFun(x, knots = knots))
    expect_equivalent(bSpline(x, degree = 2L), bsFun(x, degree = 2L))
    expect_equivalent(bSpline(x, intercept = TRUE), bsFun(x, intercept = TRUE))
    expect_equivalent(bSpline(x, knots = knots, intercept = TRUE),
                      bsFun(x, knots = knots, intercept = TRUE))
})


test_that("outputs of piecewise constant bases", {
    x <- c(NA, seq.int(0, 0.5, 0.1), NA, seq.int(0.6, 1, 0.1), NA)
    knots <- c(0.25, 0.5, 0.75)
    ## for testing splines with degree zero
    bsMat0a <- bSpline(x, degree = 0, intercept = TRUE)
    bsMat0b <- bSpline(x, df = 5, degree = 0)
    bsMat0c <- bSpline(x, df = 5, degree = 0, intercept = TRUE)
    bsMat0d <- bSpline(x, knots = knots, degree = 0)
    bsMat0e <- bSpline(x, knots = knots, degree = 0, intercept = TRUE)
    expect_output(str(bsMat0a), "matrix [1:14, 1]", fixed = TRUE)
    expect_equal(sum(is.na(bsMat0b)), 12L) # keep NA's as is
    expect_output(str(bsMat0b), "matrix [1:14, 1:4]", fixed = TRUE)
    expect_output(str(bsMat0c), "matrix [1:14, 1:5]", fixed = TRUE)
    expect_output(str(bsMat0d), "matrix [1:14, 1:3]", fixed = TRUE)
    expect_output(str(bsMat0e), "matrix [1:14, 1:4]", fixed = TRUE)
    expect_output(str(bSpline(x, df = 10, knots = knots, degree = 0L)),
                  "matrix [1:14, 1:3]", fixed = TRUE)
    expect_output(str(bSpline(x, df = 10, knots = knots,
                              degree = 0, intercept = TRUE)),
                  "matrix [1:14, 1:4]", fixed = TRUE)
    expect_error(bSpline(x, degree = 0),
                 "'intercept' has to be 'TRUE'", fixed = TRUE)
    expect_warning(bSpline(x, df = 1, knots = 0.5, degree = 0),
                   "'df' specified was not appropriate.", fixed = TRUE)
    expect_warning(bSpline(x, df = 3, knots = 0.5, degree = 0),
                   "'df' specified was not appropriate.", fixed = TRUE)
    expect_warning(bSpline(c(x, 10), knots = knots, degree = 0,
                           Boundary.knots = c(0, 1)),
                   "beyond boundary knots", fixed = TRUE)
})



## true close form formula given the all knots and degree
test_that("two internal knots, degree 0", {
    ## test with two internal knots
    x <- seq.int(0, 5, 0.1)
    knots <- c(1, 3)
    b0_1 <- function(x) as.numeric(x >= 0 & x < 1)
    b0_2 <- function(x) as.numeric(x >= 1 & x < 3)
    b0_3 <- function(x) as.numeric(x >= 3 & x < 5)
    expect_equivalent(bSpline(x, knots = knots, degree = 0L, intercept = TRUE),
                      cbind(b0_1(x), b0_2(x), b0_3(x)))
})
