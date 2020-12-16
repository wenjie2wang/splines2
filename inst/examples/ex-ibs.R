library(splines2)

x <- seq.int(0, 1, 0.01)
knots <- c(0.2, 0.4, 0.7, 0.9)
ibsMat <- ibs(x, knots = knots, degree = 1, intercept = TRUE)

## the B-spline bases integrated by function bSpline (same arguments)
bsMat0 <- bSpline(x, knots = knots, degree = 1, intercept = TRUE)

## or by the deriv method
bsMat <- deriv(ibsMat)
stopifnot(all.equal(bsMat0, bsMat, check.attributes = FALSE))

## plot B-spline basis with their corresponding integrals
op <- par(mfrow = c(1, 2))
matplot(x, bsMat, type = "l", ylab = "B-spline basis")
abline(v = knots, lty = 2, col = "gray")
matplot(x, ibsMat, type = "l", ylab = "Integral of B-spline basis")
abline(v = knots, lty = 2, col = "gray")

## reset to previous plotting settings
par(op)
