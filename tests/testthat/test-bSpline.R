context("Testing bSpline")

x <- c(NA, seq.int(0, 0.5, 0.1), NA, seq.int(0.6, 1, 0.1), NA)
knots <- c(0.25, 0.5, 0.75)

## for ease of testing
bsFun <- function(x, df, knots, degree, intercept, Boundary.knots)
{
    funCall <- match.call()
    funCall[[1L]] <- quote(bs)
    bsMat <- eval(funCall)
    class(bsMat) <- "matrix"
    bsMat
}

test_that("NA's are returned as is.", {
    expect_equivalent(bSpline(x), bsFun(x))
    expect_equivalent(bSpline(x, df = 5), bsFun(x, df = 5))
    expect_equivalent(bSpline(x, knots = knots), bsFun(x, knots = knots))
    expect_equivalent(bSpline(x, degree = 2L), bsFun(x, degree = 2L))
    expect_equivalent(bSpline(x, intercept = TRUE), bsFun(x, intercept = TRUE))
    expect_equivalent(bSpline(x, knots = knots, intercept = TRUE),
                      bsFun(x, knots = knots, intercept = TRUE))
})


## for testing splines with degree zero
bsMat0a <- bSpline(x, degree = 0, intercept = TRUE)
bsMat0b <- bSpline(x, df = 5, degree = 0)
bsMat0c <- bSpline(x, df = 5, degree = 0, intercept = TRUE)
bsMat0d <- bSpline(x, knots = knots, degree = 0)
bsMat0e <- bSpline(x, knots = knots, degree = 0, intercept = TRUE)

test_that("piece-wise constant splines", {
    expect_output(str(bsMat0a), "matrix [1:14, 1]", fixed = TRUE)
    expect_equal(sum(is.na(bsMat0b)), 12L)
    expect_output(str(bsMat0b), "matrix [1:14, 1:4]", fixed = TRUE)
    expect_output(str(bsMat0c), "matrix [1:14, 1:5]", fixed = TRUE)
    expect_output(str(bsMat0d), "matrix [1:14, 1:3]", fixed = TRUE)
    expect_output(str(bsMat0e), "matrix [1:14, 1:4]", fixed = TRUE)
    expect_output(str(bSpline(x, df = 10, knots = knots, degree = 0L)),
                  "matrix [1:14, 1:3]", fixed = TRUE)
    expect_output(str(bSpline(x, df = 10, knots = knots,
                              degree = 0, intercept = TRUE)),
                  "matrix [1:14, 1:4]", fixed = TRUE)
    expect_error(bSpline(x, degree = 0),
                 "'intercept' has to be 'TRUE'", fixed = TRUE)
    expect_warning(bSpline(x, df = 1, knots = 0.5, degree = 0),
                   "'df' specified was not appropriate.", fixed = TRUE)
    expect_warning(bSpline(x, df = 3, knots = 0.5, degree = 0),
                   "'df' specified was not appropriate.", fixed = TRUE)
    expect_warning(bSpline(c(x, 10), knots = knots, degree = 0,
                           Boundary.knots = c(0, 1)),
                   "beyond boundary knots", fixed = TRUE)
})
