library(splines2)

x <- seq.int(0, 1, 0.01)
knots <- c(0.3, 0.5, 0.6)

## quadratic B-splines
bsMat2 <- bSpline(x, knots = knots, degree = 2, intercept = TRUE)

## cubic B-splines
bsMat3 <- update(bsMat2, degree = 3)
