context("Testing cSpline")


test_that("check output dimension", {
    x <- seq.int(0, 10, 0.5)
    knots <- c(3, NA, 5, 7, NA)
    expect_equal(isNumMatrix(cSpline(x, df = 5), 21L, 5L,
                             warn_na = FALSE, error_na = FALSE), TRUE)
    expect_equal(isNumMatrix(cSpline(x, knots = knots, degree = 2),
                             21L, 5L, warn_na = FALSE, error_na = FALSE), TRUE)
    expect_equal(isNumMatrix(cSpline(x, knots = knots), 21L, 6L,
                             warn_na = FALSE, error_na = FALSE), TRUE)
    expect_equal(isNumMatrix(cSpline(x, knots = knots, intercept = TRUE),
                             21L, 7L, warn_na = FALSE, error_na = FALSE), TRUE)
    expect_error(cSpline(x, degree = 0))
    expect_equal(isNumMatrix(cSpline(x, degree = 0, intercept = TRUE), 21L, 1L,
                             warn_na = FALSE, error_na = FALSE), TRUE)
    expect_equal(isNumMatrix(cSpline(x, degree = 1), 21L, 1L,
                             warn_na = FALSE, error_na = FALSE), TRUE)
    expect_equal(isNumMatrix(cSpline(x, degree = 1, intercept = TRUE), 21L, 2L,
                             warn_na = FALSE, error_na = FALSE), TRUE)
    expect_equal(isNumMatrix(cSpline(x, degree = 2), 21L, 2L,
                             warn_na = FALSE, error_na = FALSE), TRUE)
    expect_equal(isNumMatrix(cSpline(x, degree = 2, intercept = TRUE), 21L, 3L,
                             warn_na = FALSE, error_na = FALSE), TRUE)
    expect_equal(isNumMatrix(cSpline(x), 21L, 3L,
                             warn_na = FALSE, error_na = FALSE), TRUE)
    expect_equal(isNumMatrix(cSpline(x, intercept = TRUE), 21L, 4L,
                             warn_na = FALSE, error_na = FALSE), TRUE)
})
