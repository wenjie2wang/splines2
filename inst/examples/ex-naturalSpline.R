library(splines2)

x <- seq.int(0, 1, 0.01)
knots <- c(0.3, 0.5, 0.6)

## natural spline basis
nsMat0 <- naturalSpline(x, knots = knots, intercept = TRUE)
## integrals
nsMat1 <- naturalSpline(x, knots = knots, intercept = TRUE, integral = TRUE)
## first derivatives
nsMat2 <- naturalSpline(x, knots = knots, intercept = TRUE, derivs = 1)
## second derivatives
nsMat3 <- naturalSpline(x, knots = knots, intercept = TRUE, derivs = 2)

op <- par(mfrow = c(2, 2), mar = c(2.5, 2.5, 0.2, 0.1), mgp = c(1.5, 0.5, 0))
matplot(x, nsMat0, type = "l", ylab = "basis")
matplot(x, nsMat1, type = "l", ylab = "integral")
matplot(x, nsMat2, type = "l", ylab = "1st derivative")
matplot(x, nsMat3, type = "l", ylab = "2nd derivative")
par(op) # reset to previous plotting settings

## use the deriv method
all.equal(nsMat0, deriv(nsMat1), check.attributes = FALSE)
all.equal(nsMat2, deriv(nsMat0))
all.equal(nsMat3, deriv(nsMat2))
all.equal(nsMat3, deriv(nsMat0, 2))
