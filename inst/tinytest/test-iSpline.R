## helper function
isNumMatrix <- splines2:::isNumMatrix

## check equal dimension
x <- seq.int(0, 10, 0.5)
knots <- c(3, NA, 5, 7, NA)
expect_equal(isNumMatrix(iSpline(x), 21L, 4L,
                         warn_na = FALSE, error_na = FALSE), TRUE)
expect_equal(isNumMatrix(iSpline(x, intercept = FALSE), 21L, 3L,
                         warn_na = FALSE, error_na = FALSE), TRUE)
expect_equal(iSpline(x, derivs = 1), mSpline(x, intercept = TRUE))
expect_equal(iSpline(x, df = 5, knots = knots, derivs = 1),
             mSpline(x, df = 5, knots = knots, intercept = TRUE))
expect_error(iSpline(x, knots = knots, derivs = - 1.1))
expect_warning(iSpline(x, df = 2, intercept = FALSE))
expect_equal(isNumMatrix(iSpline(x, df = 3, intercept = FALSE), 21L, 3L,
                         warn_na = FALSE, error_na = FALSE), TRUE)
expect_equal(isNumMatrix(iSpline(x, df = 5, intercept = FALSE), 21L, 5L,
                         warn_na = FALSE, error_na = FALSE), TRUE)
expect_equal(isNumMatrix(iSpline(x, knots = knots, degree = 2),
                         21L, 6L, warn_na = FALSE, error_na = FALSE),
             TRUE)
expect_equal(isNumMatrix(iSpline(x, knots = knots), 21L, 7L,
                         warn_na = FALSE, error_na = FALSE),
             TRUE)
expect_equal(isNumMatrix(iSpline(x, knots = knots, intercept = FALSE),
                         21L, 6L, warn_na = FALSE, error_na = FALSE),
             TRUE)