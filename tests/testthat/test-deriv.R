context("Testing deriv methods")

x <- c(seq.int(0, 1, 0.1), NA)
knots <- c(0.3, NA, 0.5, 0.6, NA)


test_that("test deriv methods for B-splines related", {
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
})


test_that("test deriv methods for M-splines related", {
    ## only test for scale == FALSE
    csMat <- cSpline(x, knots = knots, scale = FALSE)
    isMat <- iSpline(x, knots = knots)
    msMat <- mSpline(x, knots = knots)
    ms1Mat <- mSpline(x, knots = knots, derivs = 1)
    ms2Mat <- mSpline(x, knots = knots, derivs = 2)
    ms3Mat <- mSpline(x, knots = knots, derivs = 3)
    ms4Mat <- mSpline(x, knots = knots, derivs = 4)
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
    ## simple test for scale = TRUE
    csMat <- cSpline(x, knots = knots, degree = 4)
    expect_output(str(deriv(csMat)),
                  "matrix [1:12, 1:7]", fixed = TRUE)
    expect_equivalent(deriv(csMat, 2), deriv(deriv(csMat)))
    expect_equivalent(deriv(csMat, 3), deriv(deriv(csMat, 2)))
    expect_equivalent(deriv(csMat, 3), deriv(deriv(deriv(csMat))))
})
