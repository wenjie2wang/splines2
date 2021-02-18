library(splines2)

set.seed(123)
x <- rnorm(100)

## B-spline basis
bsMat <- bSpline(x, df = 8, degree = 3)

## extract internal knots placed based on the quantile of x
(internal_knots <- knots(bsMat))

## extract boundary knots placed based on the range of x
boundary_knots <- knots(bsMat, type = "boundary")
all.equal(boundary_knots, range(x))
