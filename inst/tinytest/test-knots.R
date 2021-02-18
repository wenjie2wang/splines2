x <- rnorm(100)

## internal knots
no_first_last <- function(x) {
    unname(x[seq.int(2, length(x) - 1)])
}
df0 <- 10
degree0 <- 3
inter_knots1 <- no_first_last(
    quantile(x, seq.int(0, 1, length.out = df0 - degree0 + 2))
)
inter_knots2 <- no_first_last(
    quantile(x, seq.int(0, 1, length.out = df0 - degree0 + 2 - 1))
)
inter_knots3 <- no_first_last(
    quantile(x, seq.int(0, 1, length.out = df0 - 1 + 2))
)

## boundary knots
boundary_knots <- range(x)

## basis matrices
bsMat <- bSpline(x, df = df0)
ibsMat <- ibs(x, df = df0)
dbsMat <- dbs(x, df = df0)
msMat <- mSpline(x, df = df0)
isMat <- iSpline(x, df = df0)
csMat <- cSpline(x, df = df0)
bpMat <- bernsteinPoly(x, df = df0)
nsMat <- naturalSpline(x, df = df0)

## internal knots
expect_equivalent(knots(bsMat), inter_knots1)
expect_equivalent(knots(ibsMat), inter_knots1)
expect_equivalent(knots(dbsMat), inter_knots1)
expect_equivalent(knots(msMat), inter_knots1)
expect_equivalent(knots(isMat), inter_knots2)
expect_equivalent(knots(csMat), inter_knots2)
expect_equivalent(knots(bpMat), NULL)
expect_equivalent(knots(nsMat), inter_knots3)
