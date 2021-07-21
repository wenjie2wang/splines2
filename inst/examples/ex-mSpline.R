library(splines2)

### example given in the reference paper by Ramsay (1988)
x <- seq.int(0, 1, 0.01)
knots <- c(0.3, 0.5, 0.6)
msMat <- mSpline(x, knots = knots, degree = 2, intercept = TRUE)

op <- par(mar = c(2.5, 2.5, 0.2, 0.1), mgp = c(1.5, 0.5, 0))
matplot(x, msMat, type = "l", ylab = "y")
abline(v = knots, lty = 2, col = "gray")

## derivatives of M-splines
dmsMat <- mSpline(x, knots = knots, degree = 2,
                  intercept = TRUE, derivs = 1)

## or using the deriv method
dmsMat1 <- deriv(msMat)
stopifnot(all.equal(dmsMat, dmsMat1, check.attributes = FALSE))

### periodic M-splines
x <- seq.int(0, 3, 0.01)
bknots <- c(0, 1)
pMat <- mSpline(x, knots = knots, degree = 3, intercept = TRUE,
                Boundary.knots = bknots, periodic = TRUE)
## integrals
iMat <- mSpline(x, knots = knots, degree = 3, intercept = TRUE,
                Boundary.knots = bknots, periodic = TRUE, integral = TRUE)
## first derivatives by "derivs = 1"
dMat1 <- mSpline(x, knots = knots, degree = 3, intercept = TRUE,
                 Boundary.knots = bknots, periodic = TRUE, derivs = 1)
## first derivatives by using the deriv() method
dMat2 <- deriv(pMat)

par(mfrow = c(2, 2))
matplot(x, pMat, type = "l", ylab = "Periodic Basis")
abline(v = seq.int(0, max(x)), lty = 2, col = "grey")
matplot(x, iMat, type = "l", ylab = "Integrals from 0")
abline(v = seq.int(0, max(x)), h = seq.int(0, max(x)), lty = 2, col = "grey")
matplot(x, dMat1, type = "l", ylab = "1st derivatives by 'derivs=1'")
abline(v = seq.int(0, max(x)), lty = 2, col = "grey")
matplot(x, dMat2, type = "l", ylab = "1st derivatives by 'deriv()'")
abline(v = seq.int(0, max(x)), lty = 2, col = "grey")
## reset to previous plotting settings
par(op)

### default placement of internal knots for periodic splines
default_knots <- function(x, df, intercept = FALSE,
                        Boundary.knots = range(x, na.rm = TRUE)) {
    ## get x in the cyclic interval [0, 1)
    x2 <- (x - Boundary.knots[1]) %% (Boundary.knots[2] - Boundary.knots[1])
    knots <- quantile(x2, probs = seq(0, 1, length.out = df + 2 - intercept))
    unname(knots[- c(1, length(knots))])
}

df <- 8
degree <- 3
intercept <- TRUE
internal_knots <- default_knots(x, df, intercept)
## 1. specify df
spline_basis1 =  splines2::mSpline(x, degree = degree, df = df,
                                   periodic = TRUE, intercept = intercept)
## 2. specify knots
spline_basis2 =  splines2::mSpline(x, degree = degree, knots = internal_knots,
                                   periodic = TRUE, intercept = intercept)

all.equal(internal_knots, knots(spline_basis1))
all.equal(spline_basis1, spline_basis2)
