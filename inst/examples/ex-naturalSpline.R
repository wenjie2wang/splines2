library(splines2)

x <- seq.int(0, 1, 0.01)
knots <- c(0.3, 0.5, 0.6)

## naturalSpline()
nsMat0 <- naturalSpline(x, knots = knots, intercept = TRUE)
nsMat1 <- naturalSpline(x, knots = knots, intercept = TRUE, integral = TRUE)
nsMat2 <- naturalSpline(x, knots = knots, intercept = TRUE, derivs = 1)
nsMat3 <- naturalSpline(x, knots = knots, intercept = TRUE, derivs = 2)

op <- par(mfrow = c(2, 2), mar = c(2.5, 2.5, 0.2, 0.1), mgp = c(1.5, 0.5, 0))
plot(nsMat0, ylab = "basis")
plot(nsMat1, ylab = "integral")
plot(nsMat2, ylab = "1st derivative")
plot(nsMat3, ylab = "2nd derivative")
par(op) # reset to previous plotting settings

## nsk()
nskMat <- nsk(x, knots = knots, intercept = TRUE)
plot(nskMat, ylab = "nsk()", mark_knots = "all")
abline(h = 1, col = "red", lty = 3)

## use the deriv method
all.equal(nsMat0, deriv(nsMat1), check.attributes = FALSE)
all.equal(nsMat2, deriv(nsMat0))
all.equal(nsMat3, deriv(nsMat2))
all.equal(nsMat3, deriv(nsMat0, 2))
