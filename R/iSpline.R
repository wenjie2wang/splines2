### I-spline bases
ims <- function(x, df = NULL, knots = NULL, degree = 3, intercept = FALSE,
               Boundary.knots = range(x), ...) {

    ## M-spline bases for inputs
    msOut <- ms(x = x, df = df, knots = knots, degree = degree,
               intercept = intercept, Boundary.knots = Boundary.knots)

    ## update input
    degree <- attr(msOut, "degree")
    knots <- attr(msOut, "knots")
    bKnots <- attr(msOut, "Boundary.knots")
    ord <- 1L + degree

    ## define knot sequence
    aKnots <- sort(c(rep(bKnots, ord + 1), knots))

    ## generate M-spline bases with (degree + 1)
    msOut1 <- ms(x = x, knots = knots, degree = ord,
                intercept = FALSE, Boundary.knots = bKnots)

    foo <- stepfun(x = knots, y = seq(ord, length(knots) + ord))
    j <- as.integer(foo(x))

    ## calculate I-spline bases at each x
    numer1 <- diff(aKnots, lag = ord + 1)[- 1L]
    msMat <- rep(numer1, each = length(x)) * msOut1 / (ord + 1)
    msAugMat <- cbind(j, msMat)
    imsOut <- t(apply(msAugMat, 1, function(b) {
        j <- b[1L]
        a <- b[- 1L]
        idx <- seq_along(a)
        a[tail(idx, - j)] <- 0
        a[seq_len(j)] <- rev(cumsum(rev(a[seq_len(j)])))
        a[idx < j - ord] <- 1
        a
    }))
    if (! intercept) imsOut <- imsOut[, - 1L, drop = FALSE]

    ## return
    list(msMat = msOut, imsMat = imsOut)
}
