x <- c(seq.int(0, 10, 0.5), NA)
knots <- c(1, 4, 8)

## bSpline
bsMat0 <- bSpline(x, degree = 2)
bsMat1 <- update(bsMat, degree = 3)
expect_equal(attr(bsMat1, "degree"), 3L)

## ibs
ibsMat0 <- ibs(x, degree = 3)
ibsMat1 <- update(ibsMat, knots = c(2, 4, 8))
expect_equal(knots(ibsMat1), c(2, 4, 8))

## dbs
dbsMat0 <- dbs(x, df = 6)
dbsMat1 <- update(dbsMat0, df = 8)
expect_equal(ncol(dbsMat1), 8L)

## mSpline
msMat0 <- mSpline(x, derivs = 1)
msMat1 <- update(msMat0, derivs = 2)
expect_equal(deriv(msMat0), msMat1)

msMat0 <- mSpline(x, integral = TRUE)
msMat1 <- update(msMat0, integral = FALSE)
expect_equivalent(deriv(msMat0), msMat1)

msMat0 <- mSpline(x, knots = knots, periodic = TRUE)
msMat1 <- update(msMat0, periodic = FALSE)
expect_false(attr(msMat1, "periodic"))

## iSpline
isMat0 <- iSpline(x, intercept = TRUE)
isMat1 <- update(isMat0, intercept = FALSE)
expect_equal(ncol(isMat1) + 1L, ncol(isMat0))

## cSpline
csMat0 <- cSpline(x)
csMat1 <- update(csMat0, scale = FALSE)
expect_true(attr(csMat0, "scale"))
expect_false(attr(csMat1, "scale"))

## bernsteinPoly
bpMat0 <- bernsteinPoly(x)
bpMat1 <- update(bpMat0, degree = 4)
expect_equal(ncol(bpMat0) + 1L, ncol(bpMat1))
expect_equal(bpMat0, update(bpMat0, knots = knots))
expect_equal(bpMat0, update(bpMat0))

## naturalSpline
nsMat0 <- naturalSpline(x, df = 3)
nsMat1 <- update(nsMat0, df = 4)
expect_equal(ncol(nsMat1), 4L)
