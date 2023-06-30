library(splines2)

set.seed(1)
x <- runif(100)
knots <- c(0.3, 0.5, 0.6) # internal knots

## cubic B-splines
bsMat <- bSpline(x, knots = knots, degree = 3, intercept = TRUE)
ibsMat <- update(bsMat, integral = TRUE) # the integrals
d1Mat <- deriv(bsMat)                    # the 1st derivaitves
d2Mat <- deriv(bsMat, 2)                 # the 2nd derivaitves

op <- par(mfrow = c(2, 2), mar = c(2.5, 2.5, 0.2, 0.1), mgp = c(1.5, 0.5, 0))
plot(bsMat, ylab = "Cubic B-splines")
plot(ibsMat, ylab = "The integrals")
plot(d1Mat, ylab = "The 1st derivatives")
plot(d2Mat, ylab = "The 2nd derivatives")

## evaluate at new values
predict(bsMat, c(0.125, 0.801))

## periodic B-splines
px <- seq(0, 3, 0.01)
pbsMat <- bSpline(px, knots = knots, Boundary.knots = c(0, 1),
                  intercept = TRUE, periodic = TRUE)
ipMat <- update(pbsMat, integral = TRUE)
dpMat <- deriv(pbsMat)
dp2Mat <- deriv(pbsMat, 2)

plot(pbsMat, ylab = "Periodic B-splines", mark_knots = "b")
plot(ipMat, ylab = "The integrals", mark_knots = "b")
plot(dpMat, ylab = "The 1st derivatives", mark_knots = "b")
plot(dp2Mat, ylab = "The 2nd derivatives", mark_knots = "b")
par(op) # reset to previous plotting settings
