library(splines2)

x <- c(seq.int(0, 1, 0.1), NA)  # NA's will be kept.
knots <- c(0.3, 0.5, 0.6)

## helper function
stopifnot_equivalent <- function(...) {
    stopifnot(all.equal(..., check.attributes = FALSE))
}

## integal of B-splines and the corresponding B-splines integrated
ibsMat <- ibs(x, knots = knots)
bsMat <- bSpline(x, knots = knots)

## the first derivative
d1Mat <- deriv(ibsMat)
stopifnot_equivalent(bsMat, d1Mat)

## the second derivative
d2Mat1 <- deriv(bsMat)
d2Mat2 <- deriv(ibsMat, derivs = 2L)
stopifnot_equivalent(d2Mat1, d2Mat2)

## nested calls are supported
d2Mat3 <- deriv(deriv(ibsMat))
stopifnot_equivalent(d2Mat2, d2Mat3)

## C-splines, I-splines, M-splines and the derivatives
csMat <- cSpline(x, knots = knots, intercept = TRUE, scale = FALSE)
isMat <- iSpline(x, knots = knots, intercept = TRUE)
stopifnot_equivalent(isMat, deriv(csMat))

msMat <- mSpline(x, knots = knots, intercept = TRUE)
stopifnot_equivalent(msMat, deriv(isMat))
stopifnot_equivalent(msMat, deriv(csMat, 2))
stopifnot_equivalent(msMat, deriv(deriv(csMat)))

dmsMat <- mSpline(x, knots = knots, intercept = TRUE, derivs = 1)
stopifnot_equivalent(dmsMat, deriv(msMat))
stopifnot_equivalent(dmsMat, deriv(isMat, 2))
stopifnot_equivalent(dmsMat, deriv(deriv(isMat)))
stopifnot_equivalent(dmsMat, deriv(csMat, 3))
stopifnot_equivalent(dmsMat, deriv(deriv(deriv(csMat))))
