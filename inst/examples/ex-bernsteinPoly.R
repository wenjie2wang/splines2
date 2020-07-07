library(splines2)

x1 <- seq.int(0, 1, 0.01)
x2 <- seq.int(- 2, 2, 0.01)

## Bernstein polynomial basis matrix over [0, 1]
bMat1 <- bernsteinPoly(x1, degree = 4, intercept = TRUE)

## generalized Bernstein polynomials basis over [- 1, 2]
bMat2 <- bernsteinPoly(x2, degree = 4, intercept = TRUE)

par(mfrow = c(1, 2), mar = c(2.5, 2.5, 0.2, 0.1), mgp = c(1.5, 0.5, 0))
matplot(x1, bMat1, type = "l", ylab = "y")
matplot(x2, bMat2, type = "l", ylab = "y")

## the first and second derivative matrix
d1Mat1 <- bernsteinPoly(x1, degree = 4, derivs = 1, intercept = TRUE)
d2Mat1 <- bernsteinPoly(x1, degree = 4, derivs = 2, intercept = TRUE)
d1Mat2 <- bernsteinPoly(x2, degree = 4, derivs = 1, intercept = TRUE)
d2Mat2 <- bernsteinPoly(x2, degree = 4, derivs = 2, intercept = TRUE)

par(mfrow = c(2, 2))
matplot(x1, d1Mat1, type = "l", ylab = "y")
matplot(x2, d1Mat2, type = "l", ylab = "y")
matplot(x1, d2Mat1, type = "l", ylab = "y")
matplot(x2, d2Mat2, type = "l", ylab = "y")

## or use the deriv method
all.equal(d1Mat1, deriv(bMat1))
all.equal(d2Mat1, deriv(bMat1, 2))

## the integrals
iMat1 <- bernsteinPoly(x1, degree = 4, integral = TRUE, intercept = TRUE)
iMat2 <- bernsteinPoly(x2, degree = 4, integral = TRUE, intercept = TRUE)
all.equal(deriv(iMat1), bMat1, check.attributes = FALSE)
all.equal(deriv(iMat2), bMat2, check.attributes = FALSE)
