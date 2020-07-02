## helper function
isNumMatrix <- splines2:::isNumMatrix

x <- seq.int(0, 10, 0.5)
knots <- c(3, 5, 7)
expect_equal(isNumMatrix(cSpline(x, df = 5), 21L, 5L,
                         warn_na = FALSE, error_na = FALSE), TRUE)
expect_equal(isNumMatrix(cSpline(x, knots = knots, degree = 2),
                         21L, 6L, warn_na = FALSE, error_na = FALSE), TRUE)
expect_equal(isNumMatrix(cSpline(x, knots = knots), 21L, 7L,
                         warn_na = FALSE, error_na = FALSE), TRUE)
expect_equal(isNumMatrix(cSpline(x, knots = knots, intercept = FALSE),
                         21L, 6L, warn_na = FALSE, error_na = FALSE), TRUE)
expect_error(cSpline(x, degree = 0, intercept = FALSE))
expect_equal(isNumMatrix(cSpline(x, degree = 0, intercept = TRUE), 21L, 1L,
                         warn_na = FALSE, error_na = FALSE), TRUE)
expect_equal(isNumMatrix(cSpline(x, degree = 1), 21L, 2L,
                         warn_na = FALSE, error_na = FALSE), TRUE)
expect_equal(isNumMatrix(cSpline(x, degree = 1, intercept = FALSE), 21L, 1L,
                         warn_na = FALSE, error_na = FALSE), TRUE)
expect_equal(isNumMatrix(cSpline(x, degree = 2), 21L, 3L,
                         warn_na = FALSE, error_na = FALSE), TRUE)
expect_equal(isNumMatrix(cSpline(x, degree = 2, intercept = FALSE), 21L, 2L,
                         warn_na = FALSE, error_na = FALSE), TRUE)
expect_equal(isNumMatrix(cSpline(x), 21L, 4L,
                         warn_na = FALSE, error_na = FALSE), TRUE)
expect_equal(isNumMatrix(cSpline(x, intercept = FALSE), 21L, 3L,
                         warn_na = FALSE, error_na = FALSE), TRUE)
