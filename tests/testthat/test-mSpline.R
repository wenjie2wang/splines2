context("Testing mSpline")


test_that("check outputs", {
    x <- c(NA, seq.int(0, 0.5, 0.1), NA, seq.int(0.6, 1, 0.1), NA)
    knots <- c(0.25, 0.5, 0.75)
    ## for testing splines with degree zero
    msMat0a <- mSpline(x, degree = 0, intercept = TRUE)
    msMat0b <- mSpline(x, knots = knots, degree = 1)
    msMat0c <- mSpline(x, knots = knots, degree = 1, intercept = TRUE)
    msMat0d <- mSpline(x, knots = knots, degree = 2)
    msMat0e <- mSpline(x, knots = knots, degree = 2, intercept = TRUE)
    expect_output(str(msMat0a), "matrix [1:14, 1]", fixed = TRUE)
    expect_equal(sum(is.na(msMat0b)), 12L) # keep NA's as is
    expect_output(str(msMat0b), "matrix [1:14, 1:4]", fixed = TRUE)
    expect_output(str(msMat0c), "matrix [1:14, 1:5]", fixed = TRUE)
    expect_output(str(msMat0d), "matrix [1:14, 1:5]", fixed = TRUE)
    expect_output(str(msMat0e), "matrix [1:14, 1:6]", fixed = TRUE)
    expect_error(mSpline(x, degree = 0),
                 "'intercept' has to be 'TRUE'", fixed = TRUE)
    expect_warning(mSpline(x, df = 1, knots = 0.5, degree = 0),
                   "'df' specified was not appropriate.", fixed = TRUE)
    expect_warning(mSpline(x, df = 3, knots = 0.5, degree = 0),
                   "'df' specified was not appropriate.", fixed = TRUE)
    expect_warning(mSpline(c(x, 10), knots = knots, degree = 0,
                           Boundary.knots = c(0, 1)),
                   "beyond boundary knots", fixed = TRUE)
})


## true close form formula given the all knots and degree
test_that("transformation of constant bases", {
    x <- seq.int(0, 7, 0.1)
    knots <- c(1, 3)
    m0_1 <- function(x) as.numeric(x < 1)
    m0_2 <- function(x) as.numeric(x >= 1 & x < 3) * 0.5
    m0_3 <- function(x) as.numeric(x >= 3 & x < 7) * 0.25
    expect_equivalent(mSpline(x, knots = knots, degree = 0L, intercept = TRUE),
                      cbind(m0_1(x), m0_2(x), m0_3(x)))
})


test_that("transformation of linear bases", {
    x <- seq.int(0, 3, 0.1)
    knots <- c(1, 2)
    ind01 <- function(x) as.numeric(x >= 0 & x < 1)
    ind12 <- function(x) as.numeric(x >= 1 & x < 2)
    ind23 <- function(x) as.numeric(x >= 2 & x <= 3)
    m0_1 <- function(x) ind01(x) * (2 - 2 * x)
    m0_2 <- function(x) ind01(x) * x + ind12(x) * (2 - x)
    m0_3 <- function(x) ind12(x) * (x - 1) + ind23(x) * (3 - x)
    m0_4 <- function(x) ind23(x) * (2 * x - 4)
    expect_equivalent(mSpline(x, knots = knots, degree = 1L, intercept = TRUE),
                      cbind(m0_1(x), m0_2(x), m0_3(x), m0_4(x)))
})
