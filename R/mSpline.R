##' M-Spline Basis for Polynomial Splines
##'
##' This function generates the M-spline basis matrix for a polynomial spline.
##'
##' It is an implementation of the close form M-spline basis based on relation
##' between M-spline basis and B-spline basis.  Internally, it calls
##' \code{\link[splines]{bs}} and generates a basis matrix for representing the
##' family of piecewise polynomials with the specified interior knots and
##' degree, evaluated at the values of \code{x}.  When "Boundary.knots" are set
##' _inside_ \code{range(x)}, \code{\link[splines]{bs}} uses a "pivot" inside
##' the respective boundary knot which is important for derivative evaluation.
##'
##' @param x The predictor variable.  Missing values are allowed.
##' @param df Degrees of freedom.  One can specify \code{df} rather than
##' \code{knots}, then the function chooses "df - degree"
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
##' @param Boundary.knots Boundary points at which to anchor the M-spline basis
##' (default the range of the non-\code{NA} data).  If both \code{knots} and
##' \code{Boundary.knots} are supplied, the basis parameters do not depend on
##' \code{x}. Data can extend beyond \code{Boundary.knots}.
##' @param ... Optional arguments for future usage.
##'
##' @return A matrix of dimension \code{length(x)} by
##' \code{df = degree + length(knots)} (plus one if intercept is included).
##' Attributes that correspond to the arguments specified are returned
##' for usage for \code{\link{predict.mSpline}}.
##'
##' @references
##' Ramsay, J. O. (1988). Monotone regression splines in action.
##' _Statistical science_, 425-441.
##'
##' @examples
##' x <- seq(0, 1, by = .01)
##' knots <- c(0.3, 0.5, 0.6)
##' mMat <- mSpline(x, knots = knots, degree = 2, intercept = TRUE)
##' matplot(x, mMat, type = "l", ylab = "M-spline bases")
##' abline(v = knots, lty = 2, col = "gray")
##'
##' @seealso
##' \code{\link{predict.mSpline}} for evaluation at given values;
##' \code{\link{iSpline}} for constructing I-spline basis.
##'
##' @importFrom splines bs
##' @export
mSpline <- function(x, df = NULL, knots = NULL, degree = 3, intercept = FALSE,
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
    class(msOut) <- c("mSpline", "basis", "matrix")
    msOut
}
