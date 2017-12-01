context("Testing bSpline")


test_that("call splines::bs", {
    x <- c(NA, seq.int(0, 0.5, 0.1), NA, seq.int(0.6, 1, 0.1), NA)
    knots <- c(0.25, NA, 0.5, 0.75, NA)
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
    knots <- c(0.25, NA, 0.5, 0.75, NA)
    ## for testing splines with degree zero
    bsMat0a <- bSpline(x, degree = 0, intercept = TRUE)
    bsMat0b <- bSpline(x, df = 5, degree = 0)
    bsMat0c <- bSpline(x, df = 5, degree = 0, intercept = TRUE)
    bsMat0d <- bSpline(x, knots = knots, degree = 0)
    bsMat0e <- bSpline(x, knots = knots, degree = 0, intercept = TRUE)
    expect_equal(isNumMatrix(bsMat0a, 14L, 1L,
                             warn_na = FALSE, error_na = FALSE), TRUE)
    expect_equal(sum(is.na(bsMat0b)), 12L) # keep NA's as is
    expect_equal(isNumMatrix(bsMat0b, 14L, 4L,
                             warn_na = FALSE, error_na = FALSE), TRUE)
    expect_equal(isNumMatrix(bsMat0c, 14L, 5L,
                             warn_na = FALSE, error_na = FALSE), TRUE)
    expect_equal(isNumMatrix(bsMat0d, 14L, 3L,
                             warn_na = FALSE, error_na = FALSE), TRUE)
    expect_equal(isNumMatrix(bsMat0e, 14L, 4L,
                             warn_na = FALSE, error_na = FALSE), TRUE)
    expect_equal(isNumMatrix(bSpline(x, df = 10, knots = knots, degree = 0L),
                             14L, 3L, warn_na = FALSE, error_na = FALSE), TRUE)
    expect_equal(isNumMatrix(bSpline(x, df = 10, knots = knots,
                                     degree = 0, intercept = TRUE),
                             14L, 4L, warn_na = FALSE, error_na = FALSE), TRUE)
    expect_error(bSpline(x, degree = 0),
                 "'intercept' has to be 'TRUE'", fixed = TRUE)
    expect_warning(bSpline(x, df = 1, knots = 0.5, degree = 0),
                   "'df' specified was not appropriate.", fixed = TRUE)
    expect_warning(bSpline(x, df = 3, knots = 0.5, degree = 0),
                   "'df' specified was not appropriate.", fixed = TRUE)
    expect_warning(bSpline(x, knots = c(- 1, 0.5), degree = 0),
                   "internal knots placed inside", fixed = TRUE)
    expect_warning(bSpline(x, knots = 0.5, degree = 0, Boundary.knots = 0:2),
                   "the first two values", fixed = TRUE)
    expect_warning(bSpline(c(x, 10), knots = knots, degree = 0,
                           Boundary.knots = c(0, 1)),
                   "beyond boundary knots", fixed = TRUE)
})



## true close form formula given the all knots and degree
test_that("two internal knots, degree 0", {
    ## test with two internal knots
    x <- seq.int(0, 5, 0.1)
    knots <- c(1, NA, 3)
    b0_1 <- function(x) as.numeric(x >= 0 & x < 1)
    b0_2 <- function(x) as.numeric(x >= 1 & x < 3)
    b0_3 <- function(x) as.numeric(x >= 3 & x < 5)
    expect_equivalent(bSpline(x, knots = knots, degree = 0L, intercept = TRUE),
                      cbind(b0_1(x), b0_2(x), b0_3(x)))
})
