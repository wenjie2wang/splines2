## get implementations of v0.2.8 for reference
v2 <- new.env()
source("../v0.2.8.R", v2)

## helper functions
isNumMatrix <- v2$isNumMatrix

### 1. check correctness first
x <- c(seq.int(0, 1, 0.1), NA)
knots <- c(0.3, 0.5, 0.6)

## test deriv methods for B-splines related
ibsMat <- ibs(x, knots = knots, intercept = TRUE)
bsMat <- bSpline(x, knots = knots, intercept = TRUE)
dbsMat <- dbs(x, knots = knots, intercept = TRUE)
dbsMatNested <- deriv(deriv(ibsMat))
d2bsMat <- dbs(x, derivs = 2, knots = knots, intercept = TRUE)
d3bsMat <- dbs(x, derivs = 3, knots = knots, intercept = TRUE)
d4bsMat <- dbs(x, derivs = 4, knots = knots, intercept = TRUE)
## first derivative of ibs == bs
expect_equivalent(bsMat, tmp <- deriv(ibsMat))
## second derivative of ibs == dbs
expect_equivalent(dbsMat, deriv(bsMat))
expect_equivalent(dbsMat, deriv(ibsMat, 2))
expect_equivalent(dbsMat, tmp <- deriv(tmp))
## third derivative of ibs == d2bs
expect_equivalent(d2bsMat, deriv(dbsMat))
expect_equivalent(d2bsMat, deriv(bsMat, 2))
expect_equivalent(d2bsMat, deriv(ibsMat, 3))
expect_equivalent(d2bsMat, tmp <- deriv(tmp))
## forth derivative of ibs == d3bs
expect_equivalent(d3bsMat, deriv(d2bsMat))
expect_equivalent(d3bsMat, deriv(dbsMat, 2))
expect_equivalent(d3bsMat, deriv(bsMat, 3))
expect_equivalent(d3bsMat, deriv(ibsMat, 4))
expect_equivalent(d3bsMat, tmp <- deriv(tmp))
## fifth derivative of ibs == d4bs
expect_equivalent(d4bsMat, deriv(d3bsMat))
expect_equivalent(d4bsMat, deriv(d2bsMat, 2))
expect_equivalent(d4bsMat, deriv(dbsMat, 3))
expect_equivalent(d4bsMat, deriv(bsMat, 4))
expect_equivalent(d4bsMat, deriv(ibsMat, 5))
expect_equivalent(d4bsMat, deriv(tmp))

## test deriv methods for M-splines related
## if scale == FALSE
csMat <- cSpline(x, knots = knots, scale = FALSE, intercept = TRUE)
isMat <- iSpline(x, knots = knots, intercept = TRUE)
msMat <- mSpline(x, knots = knots, intercept = TRUE)
ms1Mat <- mSpline(x, knots = knots, derivs = 1, intercept = TRUE)
ms2Mat <- mSpline(x, knots = knots, derivs = 2, intercept = TRUE)
ms3Mat <- mSpline(x, knots = knots, derivs = 3, intercept = TRUE)
ms4Mat <- mSpline(x, knots = knots, derivs = 4, intercept = TRUE)

## first derivative of csMat == isMat
expect_equivalent(isMat, tmp <- deriv(csMat))

## second derivative of csMat == msMat
expect_equivalent(msMat, deriv(isMat))
expect_equivalent(msMat, deriv(csMat, 2))
expect_equivalent(msMat, tmp <- deriv(tmp))

## third derivative of csMat == ms1Mat
expect_equivalent(ms1Mat, deriv(msMat))
expect_equivalent(ms1Mat, deriv(isMat, 2))
expect_equivalent(ms1Mat, deriv(csMat, 3))
expect_equivalent(ms1Mat, tmp <- deriv(tmp))

## forth derivative of csMat == ms2Mat
expect_equivalent(ms2Mat, deriv(ms1Mat))
expect_equivalent(ms2Mat, deriv(msMat, 2))
expect_equivalent(ms2Mat, deriv(isMat, 3))
expect_equivalent(ms2Mat, deriv(csMat, 4))
expect_equivalent(ms2Mat, tmp <- deriv(tmp))

## fifth derivative of csMat == ms3Mat
expect_equivalent(ms3Mat, deriv(ms2Mat))
expect_equivalent(ms3Mat, deriv(ms1Mat, 2))
expect_equivalent(ms3Mat, deriv(msMat, 3))
expect_equivalent(ms3Mat, deriv(isMat, 4))
expect_equivalent(ms3Mat, deriv(csMat, 5))
expect_equivalent(ms3Mat, tmp <- deriv(tmp))

## sixth derivative of csMat == ms4Mat
expect_equivalent(ms4Mat, deriv(ms3Mat))
expect_equivalent(ms4Mat, deriv(ms2Mat, 2))
expect_equivalent(ms4Mat, deriv(ms1Mat, 3))
expect_equivalent(ms4Mat, deriv(msMat, 4))
expect_equivalent(ms4Mat, deriv(isMat, 5))
expect_equivalent(ms4Mat, deriv(csMat, 6))
expect_equivalent(ms4Mat, tmp <- deriv(tmp))

## if scale == TRUE
csMat <- cSpline(x, knots = knots, degree = 4, scale = TRUE,
                 intercept = TRUE)
expect_true(isNumMatrix(deriv(csMat), length(x), 8L))
expect_equivalent(deriv(csMat, 2), deriv(deriv(csMat)))
expect_equivalent(deriv(csMat, 3), deriv(deriv(csMat, 2)))
expect_equivalent(deriv(csMat, 3), deriv(deriv(csMat), 2))
expect_equivalent(deriv(csMat, 3), deriv(deriv(deriv(csMat))))
expect_equivalent(deriv(csMat, 4), deriv(deriv(deriv(deriv(csMat)))))


### 2. check designed features with expectation
expect_error(deriv(ibsMat, 0), "derivs")
expect_error(deriv(bsMat, 0), "derivs")
## expect_error(deriv(csMat, 0), "derivs")
expect_error(deriv(isMat, 0), "derivs")
## expect_error(deriv(msMat, 0), "derivs")
## expect_error(deriv(dbsMat, 0), "derivs")
expect_error(deriv(ibsMat, - 1), "derivs")
expect_error(deriv(bsMat, - 1), "derivs")
expect_error(deriv(csMat, - 1), "derivs")
expect_error(deriv(isMat, - 1), "derivs")
expect_error(deriv(msMat, - 1), "derivs")
expect_error(deriv(dbsMat, - 1), "derivs")
