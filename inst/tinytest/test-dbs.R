## get implementations of v0.2.8 for reference
v2 <- new.env()
source("../v0.2.8.R", v2)

## helper functions
isNumMatrix <- v2$isNumMatrix


x <- seq.int(0, 1, 0.05)
ord <- 4
aKnots <- c(rep(0, ord), rep(1, ord))
expect_equivalent(dbs(x, derivs = 1, intercept = TRUE),
                  splines::splineDesign(aKnots, x = x, derivs = 1))
expect_equivalent(dbs(x, derivs = 2, intercept = TRUE),
                  splines::splineDesign(aKnots, x = x, derivs = 2))
## different at right boundary knot
expect_equivalent(
    dbs(x, derivs = 3, intercept = TRUE)[- length(x), ],
    splines::splineDesign(aKnots, derivs = 3, x = x)[- length(x), ]
)

knots <- c(0.2, 0.4, 0.7)
aKnots <- c(rep(0, ord), na.omit(knots), rep(1, ord))
expect_equivalent(dbs(x, derivs = 1, knots = knots, intercept = TRUE),
                  splines::splineDesign(aKnots, x = x, derivs = 1))
expect_equivalent(dbs(x, derivs = 2, knots = knots, intercept = TRUE),
                  splines::splineDesign(aKnots, x = x, derivs = 2))
## again, different at right boundary knot
expect_equivalent(
    dbs(x, derivs = 3, knots = knots, intercept = TRUE)[- length(x), ],
    splines::splineDesign(aKnots, x = x, derivs = 3)[- length(x), ]
)

knots <- c(0.3, 0.6)
ord <- 5
aKnots <- c(rep(0, ord), na.omit(knots), rep(1, ord))

expect_equivalent(dbs(x, 1, knots = knots, degree = 4, intercept = TRUE),
                  splines::splineDesign(aKnots, x, ord, derivs = 1))
expect_equivalent(dbs(x, 2, knots = knots, degree = 4, intercept = TRUE),
                  splines::splineDesign(aKnots, x, ord, derivs = 2))
expect_equivalent(dbs(x, 3, knots = knots, degree = 4, intercept = TRUE),
                  splines::splineDesign(aKnots, x, ord, derivs = 3))
## again, different at right boundary knot
expect_equivalent(
    dbs(x, 4, knots = knots, degree = 4, intercept = TRUE)[- length(x), ],
    splines::splineDesign(aKnots, x, ord, derivs = 4)[- length(x), ]
)

expect_error(dbs(x, 1, df = 1, intercept = TRUE), "df")
expect_error(dbs(x, 1, df = 2, intercept = TRUE), "df")
expect_error(dbs(x, 1, df = 3, intercept = TRUE), "df")
expect_error(dbs(x, 2, df = 3, intercept = TRUE), "df")
expect_equal(isNumMatrix(dbs(x, 1, df = 1, degree = 0, intercept = TRUE),
                         21L, 1L), TRUE)
expect_equal(isNumMatrix(dbs(x, 1, df = 4), 21L, 4L), TRUE)
expect_equal(isNumMatrix(dbs(x, 1, df = 4, intercept = TRUE),
                         21L, 4L), TRUE)
expect_equal(isNumMatrix(dbs(x, 1, df = 5),
                         21L, 5L), TRUE)
expect_equal(isNumMatrix(dbs(x, 1, df = 5, intercept = TRUE),
                         21L, 5L), TRUE)
expect_equal(isNumMatrix(dbs(x, 1, df = 5, degree = 0),
                         21L, 5L), TRUE)
expect_equal(isNumMatrix(dbs(x, 1, df = 5, degree = 0, intercept = TRUE),
                         21L, 5L), TRUE)
