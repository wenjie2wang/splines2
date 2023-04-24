source("utils.R")

x <- c(seq.int(0, 10, 0.5), NA)
bsMat <- bSpline(x, df = 6)
ibsMat <- ibs(x, df = 6)
dbsMat <- dbs(x, df = 6)
msMat <- mSpline(x, df = 6)
isMat <- iSpline(x, df = 6)
csMat1 <- cSpline(x, df = 6)
csMat2 <- cSpline(x, df = 6, scale = FALSE)
bpMat <- bernsteinPoly(x, degree = 6)
nsMat <- naturalSpline(x, df = 6)

## with newx
expect_eqt(predict(bsMat, 1), bsMat[3L, , drop = FALSE])
expect_eqt(predict(ibsMat, 1), ibsMat[3L, , drop = FALSE])
expect_eqt(predict(dbsMat, 1), dbsMat[3L, , drop = FALSE])
expect_eqt(predict(msMat, 1), msMat[3L, , drop = FALSE])
expect_eqt(predict(isMat, 1), isMat[3L, , drop = FALSE])
expect_eqt(predict(csMat1, 1), csMat1[3L, , drop = FALSE])
expect_eqt(predict(csMat2, 1), csMat2[3L, , drop = FALSE])
expect_eqt(predict(bpMat, 1), bpMat[3L, , drop = FALSE])
expect_eqt(predict(nsMat, 1), nsMat[3L, , drop = FALSE])

## without newx
expect_eqt(predict(bsMat), bsMat)
expect_eqt(predict(ibsMat), ibsMat)
expect_eqt(predict(dbsMat), dbsMat)
expect_eqt(predict(msMat), msMat)
expect_eqt(predict(isMat), isMat)
expect_eqt(predict(csMat1), csMat1)
expect_eqt(predict(csMat2), csMat2)
expect_eqt(predict(bpMat), bpMat)
expect_eqt(predict(nsMat), nsMat)

## with coef
beta <- runif(ncol(bsMat))
expect_eqt(predict(bsMat, coef = beta), as.numeric(bsMat %*% beta))
expect_eqt(predict(ibsMat, coef = beta), as.numeric(ibsMat %*% beta))
expect_eqt(predict(dbsMat, coef = beta), as.numeric(dbsMat %*% beta))
expect_eqt(predict(msMat, coef = beta), as.numeric(msMat %*% beta))
expect_eqt(predict(isMat, coef = beta), as.numeric(isMat %*% beta))
expect_eqt(predict(csMat1, coef = beta), as.numeric(csMat1 %*% beta))
expect_eqt(predict(csMat2, coef = beta), as.numeric(csMat2 %*% beta))
expect_eqt(predict(bpMat, coef = beta), as.numeric(bpMat %*% beta))
expect_eqt(predict(nsMat, coef = beta), as.numeric(nsMat %*% beta))

## update without coef
expect_eqt(predict(bsMat, derivs = 1), update(bsMat, derivs = 1))
expect_eqt(predict(ibsMat, derivs = 1), update(ibsMat, derivs = 1))
expect_eqt(predict(dbsMat, derivs = 1), update(dbsMat, derivs = 1))
expect_eqt(predict(msMat, derivs = 1), update(msMat, derivs = 1))
expect_eqt(predict(isMat, derivs = 1), update(isMat, derivs = 1))
expect_eqt(predict(csMat1, derivs = 1), update(csMat1, derivs = 1))
expect_eqt(predict(csMat2, derivs = 1), update(csMat2, derivs = 1))
expect_eqt(predict(bpMat, derivs = 1), update(bpMat, derivs = 1))
expect_eqt(predict(nsMat, derivs = 1), update(nsMat, derivs = 1))

## update with coef
expect_eqt(predict(bsMat, derivs = 1, coef = beta),
           as.numeric(deriv(bsMat, derivs = 1) %*% beta))
expect_eqt(predict(ibsMat, derivs = 1, coef = beta),
           as.numeric(deriv(ibsMat, derivs = 1) %*% beta))
expect_eqt(predict(dbsMat, derivs = 1, coef = beta),
           as.numeric(deriv(dbsMat, derivs = 1) %*% beta))
expect_eqt(predict(msMat, derivs = 1, coef = beta),
           as.numeric(deriv(msMat, derivs = 1) %*% beta))
expect_eqt(predict(isMat, derivs = 1, coef = beta),
           as.numeric(deriv(isMat, derivs = 1) %*% beta))
expect_eqt(predict(csMat1, derivs = 1, coef = beta),
           as.numeric(deriv(csMat1, derivs = 1) %*% beta))
expect_eqt(predict(csMat2, derivs = 1, coef = beta),
           as.numeric(deriv(csMat2, derivs = 1) %*% beta))
expect_eqt(predict(bpMat, derivs = 1, coef = beta),
           as.numeric(deriv(bpMat, derivs = 1) %*% beta))
expect_eqt(predict(nsMat, derivs = 1, coef = beta),
           as.numeric(deriv(nsMat, derivs = 1) %*% beta))
