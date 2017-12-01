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
    msMat0f <- mSpline(0.1, knots = knots, degree = 2, intercept = TRUE,
                       Boundary.knots = c(0, 1), derivs = 1)
    msMat0g <- mSpline(0.1, knots = knots, degree = 2, intercept = TRUE,
                       Boundary.knots = c(0, 1), derivs = 2)
    msMat0h <- mSpline(0.1, knots = knots, degree = 2,
                       Boundary.knots = c(0, 1), derivs = 3)
    expect_equal(isNumMatrix(msMat0a, 14L, 1L,
                             warn_na = FALSE, error_na = FALSE), TRUE)
    expect_equal(sum(is.na(msMat0b)), 12L) # keep NA's as is
    expect_equal(isNumMatrix(msMat0b, 14L, 4L,
                             warn_na = FALSE, error_na = FALSE), TRUE)
    expect_equal(isNumMatrix(msMat0c, 14L, 5L,
                             warn_na = FALSE, error_na = FALSE), TRUE)
    expect_equal(isNumMatrix(msMat0d, 14L, 5L,
                             warn_na = FALSE, error_na = FALSE), TRUE)
    expect_equal(isNumMatrix(msMat0e, 14L, 6L,
                             warn_na = FALSE, error_na = FALSE), TRUE)
    expect_equal(isNumMatrix(msMat0f, 1L, 6L,
                             warn_na = FALSE, error_na = FALSE), TRUE)
    expect_equal(isNumMatrix(msMat0g, 1L, 6L,
                             warn_na = FALSE, error_na = FALSE), TRUE)
    expect_equal(isNumMatrix(msMat0h, 1L, 5L,
                             warn_na = FALSE, error_na = FALSE), TRUE)
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
