library(splines2)

x <- seq.int(0, 1, 0.01)
knots <- c(0.3, 0.5, 0.6)

## cubic B-splines
bsMat <- bSpline(x, knots = knots, degree = 3, intercept = TRUE)

par(mar = c(2.5, 2.5, 0.2, 0.1), mgp = c(1.5, 0.5, 0))
matplot(x, bsMat, type = "l", ylab = "Cubic B-spline Bases")
abline(v = knots, lty = 2, col = "gray")

## the first derivaitves
d1Mat <- deriv(bsMat)

## the second derivaitves
d2Mat <- deriv(bsMat, 2)

## evaluate at new values
predict(bsMat, c(0.125, 0.801))
