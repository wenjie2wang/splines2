library(splines2)

x <- seq(0, 1, 0.1)
## Bernstein polynomial basis matrix
bMat <- bernsteinPoly(x, degree = 4, intercept = TRUE)

## the first and second derivative matrix
d1Mat <- bernsteinPoly(x, degree = 4, derivs = 1, intercept = TRUE)
d2Mat <- bernsteinPoly(x, degree = 4, derivs = 2, intercept = TRUE)
## or use the deriv method
all.equal(d1Mat, deriv(bMat))
all.equal(d2Mat, deriv(bMat, 2))

## the integral matrix
iMat <- bernsteinPoly(x, degree = 4, integral = TRUE, intercept = TRUE)
all.equal(deriv(iMat), bMat, check.attributes = FALSE)
