library(splines2)

x <- seq.int(0, 1, 0.01)
knots <- c(0.3, 0.5, 0.6)

### when 'scale = TRUE' (by default)
csMat <- cSpline(x, knots = knots, degree = 2)

op <- par(mar = c(2.5, 2.5, 0.2, 0.1), mgp = c(1.5, 0.5, 0))
matplot(x, csMat, type = "l", ylab = "C-spline basis")
abline(v = knots, lty = 2, col = "gray")
isMat <- deriv(csMat)
msMat <- deriv(csMat, derivs = 2)
matplot(x, isMat, type = "l", ylab = "scaled I-spline basis")
matplot(x, msMat, type = "l", ylab = "scaled M-spline basis")

## reset to previous plotting settings
par(op)

### when 'scale = FALSE'
csMat <- cSpline(x, knots = knots, degree = 2, scale = FALSE)

## the corresponding I-splines and M-splines (with same arguments)
isMat <- iSpline(x, knots = knots, degree = 2)
msMat <- mSpline(x, knots = knots, degree = 2, intercept = TRUE)

## or using deriv methods (more efficient)
isMat1 <- deriv(csMat)
msMat1 <- deriv(csMat, derivs = 2)

## equivalent
stopifnot(all.equal(isMat, isMat1, check.attributes = FALSE))
stopifnot(all.equal(msMat, msMat1, check.attributes = FALSE))
