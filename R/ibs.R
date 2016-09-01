##' Integral of B-Spline Basis for Polynomial Splines
##'
##' This function generates the integral of B-spline basis matrix
##' for a polynomial spline. The arguments are exactly the same with function
##' \code{\link[splines]{bs}} in package \code{\link{splines}}.
##'
##' It is an implementation of the close form integral of B-spline basis based
##' on recursion relation.  Internally, it calls \code{\link[splines]{bs}} and
##' generates a basis matrix for representing the family of piecewise
##' polynomials and their corresponding integrals with the specified interior
##' knots and degree, evaluated at the values of \code{x}.
##' When "Boundary.knots" are set _inside_ \code{range(x)},
##' \code{\link[splines]{bs}} uses a "pivot" inside the respective boundary
##' knot which is important for derivative evaluation.
##'
##' @param x The predictor variable.  Missing values are allowed.
##' @param df Degrees of freedom.  One can specify \code{df} rather than
##' \code{knots}, then \code{intBs} chooses "df - degree"
##' (minus one if there is an intercept) knots at suitable quantiles of \code{x}
##' (which will ignore missing values).  The default, \code{NULL}, corresponds
##' to _no_ inner knots, i.e., "degree - intercept".
##' @param knots The _internal_ breakpoints that define the spline.  The
##' default is \code{NULL}, which results in a basis for ordinary
##' polynomial regression.  Typical values are the mean or median
##' for one knot, quantiles for more knots.  See also
##' \code{Boundary.knots}.
##' @param degree Degree of the piecewise polynomial. The default value is 3
##' for cubic splines.
##' @param intercept If \code{TRUE}, an intercept is included in the basis;
##' default is \code{FALSE}.
##' @param Boundary.knots boundary points at which to anchor the B-spline basis
##' (default the range of the non-\code{NA} data).  If both \code{knots} and
##' \code{Boundary.knots} are supplied, the basis parameters do not depend on
##' \code{x}. Data can extend beyond \code{Boundary.knots}.
##'
##' @return A \code{ibs} class object that includes the origin B-spline and its
##' integral basis matrix. The dimension of each matrix is \code{length(x)} by
##' \code{df = degree + length(knots)} (plus on if intercept is included).
##'
##' @references
##'
##'
##'
##' @examples
##'
##'
##'
##' @seealso
##' \code{\link{ms}} for constructing M-spline basis;
##' \code{\link{ims}} for constructing I-spine basis.
##'
##' @importFrom splines bs
##' @export
ibs <- function(x, df = NULL, knots = NULL, degree = 3, intercept = FALSE,
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

    ## generate B-spline bases with (degree + 1)
    bsOut1 <- splines::bs(x = x, knots = knots, degree = ord,
                         intercept = FALSE, Boundary.knots = bKnots)
    numer1 <- diff(aKnots, lag = ord)
    if (! intercept) {
        bsOut1 <- bsOut1[, - 1L, drop = FALSE]
        numer1 <- numer1[- 1L]
    }
    numer2 <- apply(bsOut1, 1, function(a) rev(cumsum(rev(a))))
    ibsOut <- t(numer1 * numer2) / ord

    ## return
    out <- list(bsMat = bsOut, ibsMat = ibsOut)
    class(out) <- "ibs"
    out
}
