context("Testing iSpline")


test_that("check output dimension", {
    x <- seq.int(0, 10, 0.5)
    knots <- c(3, NA, 5, 7, NA)
    expect_output(str(iSpline(x)),
                  "matrix [1:21, 1:3]", fixed = TRUE)
    expect_output(str(iSpline(x, intercept = TRUE)),
                  "matrix [1:21, 1:4]", fixed = TRUE)
    expect_equal(iSpline(x, derivs = 1), mSpline(x))
    expect_equal(iSpline(x, df = 5, knots = knots, derivs = 1),
                 mSpline(x, df = 5, knots = knots))
    expect_error(iSpline(x, knots = knots, derivs = - 1.1))
    expect_warning(iSpline(x, df = 2))
    expect_output(str(iSpline(x, df = 3)),
                  "matrix [1:21, 1:3]", fixed = TRUE)
    expect_output(str(iSpline(x, df = 5)),
                  "matrix [1:21, 1:5]", fixed = TRUE)
    expect_output(str(iSpline(x, knots = knots, degree = 2)),
                  "matrix [1:21, 1:5]", fixed = TRUE)
    expect_output(str(iSpline(x, knots = knots)),
                  "matrix [1:21, 1:6]", fixed = TRUE)
    expect_output(str(iSpline(x, knots = knots, intercept = TRUE)),
                  "matrix [1:21, 1:7]", fixed = TRUE)
    expect_output(str(iSpline(x, knots = knots, intercept = TRUE)),
                  "matrix [1:21, 1:7]", fixed = TRUE)
})
