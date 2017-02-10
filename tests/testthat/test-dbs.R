context("Testing dbs with splines::splineDesign")


test_that("cubic B-splines without internal knots", {
    x <- seq.int(0, 1, 0.05)
    ord <- 4
    aKnots <- c(rep(0, ord), rep(1, ord))
    expect_equivalent(dbs(x, derivs = 1, intercept = TRUE),
                      splines::splineDesign(aKnots, x = x, derivs = 1))
    expect_equivalent(dbs(x, derivs = 2, intercept = TRUE),
                      splines::splineDesign(aKnots, x = x, derivs = 2))
    expect_equivalent(dbs(x, derivs = 3, intercept = TRUE),
                      splines::splineDesign(aKnots, x = x, derivs = 3))
})


test_that("cubic B-splines with three internal knots", {
    x <- seq.int(0, 1, 0.05)
    knots <- c(0.2, NA, 0.4, 0.7, NA)
    ord <- 4
    aKnots <- c(rep(0, ord), na.omit(knots), rep(1, ord))
    expect_equivalent(dbs(x, derivs = 1, knots = knots, intercept = TRUE),
                      splines::splineDesign(aKnots, x = x, derivs = 1))
    expect_equivalent(dbs(x, derivs = 2, knots = knots, intercept = TRUE),
                      splines::splineDesign(aKnots, x = x, derivs = 2))
    expect_equivalent(dbs(x, derivs = 3, knots = knots, intercept = TRUE),
                      splines::splineDesign(aKnots, x = x, derivs = 3))
})


test_that("quad B-splines with two internal knots", {
    x <- seq.int(0, 1, 0.05)
    knots <- c(0.3, NA, 0.6)
    ord <- 5
    aKnots <- c(rep(0, ord), na.omit(knots), rep(1, ord))
    expect_equivalent(dbs(x, 1, knots = knots, degree = 4, intercept = TRUE),
                      splines::splineDesign(aKnots, x, ord, derivs = 1))
    expect_equivalent(dbs(x, 2, knots = knots, degree = 4, intercept = TRUE),
                      splines::splineDesign(aKnots, x, ord, derivs = 2))
    expect_equivalent(dbs(x, 3, knots = knots, degree = 4, intercept = TRUE),
                      splines::splineDesign(aKnots, x, ord, derivs = 3))
    expect_equivalent(dbs(x, 4, knots = knots, degree = 4, intercept = TRUE),
                      splines::splineDesign(aKnots, x, ord, derivs = 4))
})


test_that("B-splines derivatives with df given", {
    x <- seq.int(0, 1, 0.05)
    expect_warning(dbs(x, 1, df = 0, intercept = TRUE))
    expect_warning(dbs(x, 1, df = 1, intercept = TRUE))
    expect_warning(dbs(x, 1, df = 2, intercept = TRUE))
    expect_warning(dbs(x, 1, df = 3, intercept = TRUE))
    expect_warning(dbs(x, 2, df = 3, intercept = TRUE))
    expect_output(str(dbs(x, 1, df = 1, degree = 0, intercept = TRUE)),
                  "matrix [1:21, 1]", fixed = TRUE)
    expect_output(str(dbs(x, 1, df = 4)),
                  "matrix [1:21, 1:4]", fixed = TRUE)
    expect_output(str(dbs(x, 1, df = 4, intercept = TRUE)),
                  "matrix [1:21, 1:4]", fixed = TRUE)
    expect_output(str(dbs(x, 1, df = 5)),
                  "matrix [1:21, 1:5]", fixed = TRUE)
    expect_output(str(dbs(x, 1, df = 5, intercept = TRUE)),
                  "matrix [1:21, 1:5]", fixed = TRUE)
    expect_output(str(dbs(x, 1, df = 5, degree = 0)),
                  "matrix [1:21, 1:5]", fixed = TRUE)
    expect_output(str(dbs(x, 1, df = 5, degree = 0, intercept = TRUE)),
                  "matrix [1:21, 1:5]", fixed = TRUE)
})
