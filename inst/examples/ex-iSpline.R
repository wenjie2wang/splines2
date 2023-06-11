library(splines2)

## an example given in Ramsay (1988)
x <- seq.int(0, 1, by = 0.01)
knots <- c(0.3, 0.5, 0.6)
isMat <- iSpline(x, knots = knots, degree = 2)

op <- par(mar = c(2.5, 2.5, 0.2, 0.1), mgp = c(1.5, 0.5, 0))
plot(isMat, ylab = "I-spline basis", mark_knots = "internal")
par(op) # reset to previous plotting settings

## the derivative of I-splines is M-spline
msMat1 <- iSpline(x, knots = knots, degree = 2, derivs = 1)
msMat2 <- mSpline(x, knots = knots, degree = 2, intercept = TRUE)
stopifnot(all.equal(msMat1, msMat2))
