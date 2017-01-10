context("Testing iSpline")


test_that("check output dimension", {
    x <- seq.int(0, 10, 0.5)
    knots <- c(3, 5, 7)
    expect_output(str(iSpline(x, df = 5)),
                  "matrix [1:21, 1:5]", fixed = TRUE)
    expect_output(str(iSpline(x, knots = knots, degree = 2)),
                  "matrix [1:21, 1:5]", fixed = TRUE)
    expect_output(str(iSpline(x, knots = knots)),
                  "matrix [1:21, 1:6]", fixed = TRUE)
    expect_output(str(iSpline(x, knots = knots, intercept = TRUE)),
                  "matrix [1:21, 1:7]", fixed = TRUE)
})
