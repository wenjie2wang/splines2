### M-spline
ms <- function(x, df = NULL, knots = NULL, degree = 3, intercept = FALSE,
              Boundary.knots = range(x), ...) {

    ## B-spline bases for inputs
    bsOut <- splines::bs(x = x, df = df, knots = knots, degree = degree,
                        intercept = intercept, Boundary.knots = Boundary.knots)

    ## update input
    degree <- attr(bsOut, "degree")
    knots <- attr(bsOut, "knots")
    bKnots <- attr(bsOut, "Boundary.knots")
    ord <- 1L + degree

    ## define knot sequence
    aKnots <- sort(c(rep(bKnots, ord), knots))

    ## transformation from B-splines to M-splines
    denom <- diff(aKnots, lag = ord)
    transCoef <- ifelse(abs(denom) < .Machine$double.eps, 0, ord / denom)
    if (! intercept) transCoef <- transCoef[- 1L]
    msOut <- rep(transCoef, each = length(x)) * bsOut
    attr(msOut, "class") <- NULL

    ## return
    msOut
}
