## test the $ method
x <- c(NA, seq.int(0, 0.5, 0.1), NA, seq.int(0.6, 1, 0.1), NA)
knots <- c(0.25, 0.5, 0.75)
b_knots <- c(0, 1)

## periodic = FALSE
res <- bSpline(x, df = 6)
plot(res, mark_knots = "none") # current default
plot(res, mark_knots = "internal")
plot(res, mark_knots = "boundary")
plot(res, mark_knots = "all")

## periodic = TRUE
res <- bSpline(c(x, 2 * x), df = 6, periodic = TRUE)
plot(res, mark_knots = "none") # current default
plot(res, mark_knots = "internal")
plot(res, mark_knots = "boundary")
plot(res, mark_knots = "all")

## test warning messages
x <- seq.int(- 1, 2, by = 0.1)
bs0 <- bsp(x, degree = 0, knots = 0.4, Boundary.knots = c(0, 1),
           intercept = TRUE)
expect_silent(plot(bs0))
