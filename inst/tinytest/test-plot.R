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
