context("Testing cSpline")


test_that("check output dimension", {
    x <- seq.int(0, 10, 0.5)
    knots <- c(3, NA, 5, 7, NA)
    expect_output(str(cSpline(x, df = 5)),
                  "matrix [1:21, 1:5]", fixed = TRUE)
    expect_output(str(cSpline(x, knots = knots, degree = 2)),
                  "matrix [1:21, 1:5]", fixed = TRUE)
    expect_output(str(cSpline(x, knots = knots)),
                  "matrix [1:21, 1:6]", fixed = TRUE)
    expect_output(str(cSpline(x, knots = knots, intercept = TRUE)),
                  "matrix [1:21, 1:7]", fixed = TRUE)
    expect_error(cSpline(x, degree = 0))
    expect_output(str(cSpline(x, degree = 0, intercept = TRUE)),
                  "matrix [1:21, 1]", fixed = TRUE)
    expect_output(str(cSpline(x, degree = 1)),
                  "matrix [1:21, 1]", fixed = TRUE)
    expect_output(str(cSpline(x, degree = 1, intercept = TRUE)),
                  "matrix [1:21, 1:2]", fixed = TRUE)
    expect_output(str(cSpline(x, degree = 2)),
                  "matrix [1:21, 1:2]", fixed = TRUE)
    expect_output(str(cSpline(x, degree = 2, intercept = TRUE)),
                  "matrix [1:21, 1:3]", fixed = TRUE)
    expect_output(str(cSpline(x)),
                  "matrix [1:21, 1:3]", fixed = TRUE)
    expect_output(str(cSpline(x, intercept = TRUE)),
                  "matrix [1:21, 1:4]", fixed = TRUE)
})
