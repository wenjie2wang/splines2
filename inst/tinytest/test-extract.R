## test the $ method
x <- c(NA, seq.int(0, 0.5, 0.1), NA, seq.int(0.6, 1, 0.1), NA)
knots <- c(0.25, 0.5, 0.75)
b_knots <- c(0, 1)

res <- bSpline(x, df = 6)
expect_equal(res$knots, knots(res))
expect_equal(res$Boundary.knots, b_knots)
expect_equal(res$degree, 3L)
expect_equal(res$derivs, 0L)
expect_false(res$integral)
expect_false(res$periodic)
expect_false(res$intercept)
