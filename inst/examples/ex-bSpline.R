library(splines2)
x <- seq.int(0, 1, 0.01)
knots <- c(0.3, 0.5, 0.6)
bsMat <- bSpline(x, knots = knots, degree = 0, intercept = TRUE)

library(graphics)
matplot(x, bsMat, type = "l", ylab = "Piecewise constant B-spline bases")
abline(v = knots, lty = 2, col = "gray")
