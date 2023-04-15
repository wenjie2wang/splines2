library(splines2)

x <- runif(100)
knots <- c(0.3, 0.5, 0.6) # internal knots

## cubic B-splines
bsMat <- bSpline(x, knots = knots, degree = 3, intercept = TRUE)
## the integrals
ibsMat <- ibs(x, knots = knots, degree = 3, intercept = TRUE)
## the first derivatives
d1Mat <- deriv(bsMat)
## the second derivaitves
d2Mat <- deriv(bsMat, 2)

op <- par(mfrow = c(2, 2))
plot(bsMat, main = "Cubic B-splines")
plot(ibsMat, main = "The integrals")
plot(d1Mat, main = "The first derivatives")
plot(d2Mat, main = "The second derivatives")
par(op) # reset to previous plotting settings

## evaluate at new values
predict(bsMat, c(0.125, 0.801))

## periodic B-splines
px <- seq(0, 3, 0.01)
pbsMat <- bSpline(px, knots = knots, Boundary.knots = c(0, 1),
                  intercept = TRUE, periodic = TRUE)
ipMat <- ibs(px, knots = knots, Boundary.knots = c(0, 1),
             intercept = TRUE, periodic = TRUE)
dpMat <- deriv(pbsMat)

op <- par(mfrow = c(1, 3))
plot(pbsMat, main = "Periodic B-splines", mark_knots = "b")
plot(ipMat, main = "The integrals", mark_knots = "b")
plot(dpMat, main = "The first derivatives", mark_knots = "b")
par(op) # reset to previous plotting settings
