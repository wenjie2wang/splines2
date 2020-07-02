library(splines2)
x <- seq.int(0, 1, 0.01)
knots <- c(0.2, 0.4, 0.7)
## the second derivative of cubic B-splines with three internal knots
dMat <- dbs(x, derivs = 2L, knots = knots, intercept = TRUE)

## compare with the results from splineDesign
ord <- attr(dMat, "degree") + 1L
bKnots <- attr(dMat, "Boundary.knots")
aKnots <- c(rep(bKnots[1L], ord), knots, rep(bKnots[2L], ord))
res <- splines::splineDesign(aKnots, x = x, derivs = 2L)
stopifnot(all.equal(res, dMat, check.attributes = FALSE))
