##' B-Spline Basis for Polynomial Splines
##'
##' Generates the B-spline basis matrix representing the family of piecewise
##' polynomials with the specified interior knots and degree, evaluated at the
##' values of \code{x}.
##'
##' It is an augmented function of \code{bs} function in \code{splines} package
##' for B-spline basis that allows piecewise constant ( left-closed and
##' right-open except on the right boundary) spline basis with zero degree.
##'
##' @param x The predictor variable.  Missing values are allowed and will be
##'     returned as they were.
##' @param df Degrees of freedom.  One can specify \code{df} rather than
##'     \code{knots}, then the function chooses "df - degree" (minus one if
##'     there is an intercept) knots at suitable quantiles of \code{x} (which
##'     will ignore missing values).  The default, \code{NULL}, corresponds to
##'     no inner knots, i.e., "degree - intercept". If \code{knots} was
##'     specified, \code{df} specified will be ignored.
##' @param knots The internal breakpoints that define the spline.  The default
##'     is \code{NULL}, which results in a basis for ordinary polynomial
##'     regression.  Typical values are the mean or median for one knot,
##'     quantiles for more knots.  See also \code{Boundary.knots}.
##' @param degree Non-negative integer degree of the piecewise polynomial. The
##'     default value is 3 for cubic splines. Zero degree is allowed for this
##'     function, which is the only difference compared with
##'     \code{\link[splines]{bs}} in package \code{splines}.
##' @param intercept If \code{TRUE}, an intercept is included in the basis;
##'     Default is \code{FALSE}.
##' @param Boundary.knots Boundary points at which to anchor the B-spline basis.
##'     By default, they are the range of the non-\code{NA} data.  If both
##'     \code{knots} and \code{Boundary.knots} are supplied, the basis
##'     parameters do not depend on \code{x}. Data can extend beyond
##'     \code{Boundary.knots}.
##' @param ... Optional arguments for future usage.
##'
##' @return A matrix of dimension \code{length(x)} by
##' \code{df = degree + length(knots)} (plus one if intercept is included).
##' Attributes that correspond to the arguments specified are returned
##' for usage of other functions in this package.
##'
##' @examples
##' library(splines2)
##' x <- seq.int(0, 1, 0.01)
##' knots <- c(0.3, 0.5, 0.6)
##' bsMat <- bSpline(x, knots = knots, degree = 0, intercept = TRUE)
##'
##' library(graphics)
##' matplot(x, bsMat, type = "l", ylab = "Piecewise constant B-spline bases")
##' abline(v = knots, lty = 2, col = "gray")
##'
##' @seealso
##' \code{\link{predict.bSpline2}} for evaluation at given (new) values;
##' \code{\link{dbs}}, \code{\link{deriv.bSpline2}} for derivatives;
##' \code{\link{ibs}} for integral of B-splines;
##' \code{\link{mSpline}} for M-splines;
##' \code{\link{iSpline}} for I-splines;
##' \code{\link{cSpline}} for C-splines.
##'
##' @export
bSpline <- function(x, df = NULL, knots = NULL, degree = 3L,
                    intercept = FALSE, Boundary.knots = NULL, ...)
{
    ## check inputs
    if ((degree <- as.integer(degree)) < 0)
        stop("'degree' must be a nonnegative integer.")
    if (is.null(df)) {
        df <- 0L
    } else {
        df <- as.integer(df)
        if (df < 0) {
            stop("'df' must be a nonnegative integer.")
        }
    }
    knots <- null2num0(knots)
    Boundary.knots <- null2num0(Boundary.knots)
    ## take care of possible NA's in `x`
    nax <- is.na(x)
    if (all(nax)) {
        stop("The 'x' cannot be all NA's!")
    }
    ## remove NA's in x
    xx <- if (nas <- any(nax)) {
              x[! nax]
          } else {
              x
          }
    ## call the engine function
    out <- rcpp_bSpline_basis(
        x = xx,
        df = df,
        degree = degree,
        internal_knots = knots,
        boundary_knots = Boundary.knots,
        complete_basis = intercept
    )
    ## keep NA's as is
    if (nas) {
        nmat <- matrix(NA_real_, length(nax), ncol(out))
        nmat[! nax, ] <- out
        saved_attr <- attributes(out)
        saved_attr$dim[1] <- length(nax)
        out <- nmat
        attributes(out) <- saved_attr
        attr(out, "x") <- x
    }
    ## add dimnames for consistency with bs returns
    name_x <- names(x)
    if (! is.null(name_x)) {
        row.names(out) <- name_x
    }
    ## add class
    class(out) <- c("matrix", "bSpline2")
    ## return
    out
}
