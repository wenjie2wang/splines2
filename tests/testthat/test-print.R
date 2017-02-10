context("Testing print methods")


test_that("print methods", {
    x <- c(NA, seq.int(0, 10, 0.5), NA)
    bsMat <- bSpline(x)
    ibsMat <- ibs(x)
    dbsMat <- dbs(x)
    msMat <- mSpline(x)
    isMat <- iSpline(x)
    csMat <- cSpline(x)
    expect_equivalent(print(bsMat), bsMat)
    expect_equivalent(print(ibsMat), ibsMat)
    expect_equivalent(print(dbsMat), dbsMat)
    expect_equivalent(print(msMat), msMat)
    expect_equivalent(print(isMat), isMat)
    expect_equivalent(print(csMat), csMat)
})
