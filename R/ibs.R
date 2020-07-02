##' Integral of B-Spline Basis
##'
##' Generates the integral of B-spline basis matrix for a polynomial spline.
##'
##' It is an implementation of the close form integral of B-spline basis based
##' on recursion relation.  Internally, it calls \code{\link{bSpline}} and
##' generates a basis matrix for representing the family of piecewise
##' polynomials and their corresponding integrals with the specified interior
##' knots and degree, evaluated at the values of \code{x}.
##'
##' @inheritParams bSpline
##'
##' @return
##' A matrix of dimension \code{length(x)} by
##' \code{df = degree + length(knots)} (plus on if intercept is included).
##' Attributes that correspond to the arguments specified are returned
##' for usage of other functions in this package.
##'
##' @references
##' De Boor, Carl. (1978). \emph{A practical guide to splines}.
##' Vol. 27. New York: Springer-Verlag.
##'
##' @examples
##' library(splines2)
##' x <- seq.int(0, 1, 0.01)
##' knots <- c(0.2, 0.4, 0.7, 0.9)
##' ibsMat <- ibs(x, knots = knots, degree = 1, intercept = TRUE)
##'
##' ## the B-spline bases integrated by function bSpline (same arguments)
##' bsMat0 <- bSpline(x, knots = knots, degree = 1, intercept = TRUE)
##' ## or by function deriv (recommended) that directly extracts the existing
##' ## result from the attribute of ibsMat and thus is much more efficient.
##' bsMat <- deriv(ibsMat)
##' stopifnot(all.equal(bsMat0, bsMat, check.attributes = FALSE)) # equivalent
##'
##' ## plot B-spline basis with their corresponding integrals
##' library(graphics)
##' par(mfrow = c(1, 2))
##' matplot(x, bsMat, type = "l", ylab = "B-spline basis")
##' abline(v = knots, lty = 2, col = "gray")
##' matplot(x, ibsMat, type = "l", ylab = "Integral of B-spline basis")
##' abline(v = knots, lty = 2, col = "gray")
##' par(mfrow = c(1, 1))
##'
##' @seealso
##' \code{\link{predict.ibs}} for evaluation at given (new) values;
##' \code{\link{deriv.ibs}} for derivative method.
##' \code{\link{bSpline}} for B-splines;
##' \code{\link{dbs}} for derivatives of B-splines;
##'
##' @export
ibs <- function(x, df = NULL, knots = NULL, degree = 3,
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
    out <- rcpp_bSpline_integral(
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
    class(out) <- c("matrix", "ibs")
    ## return
    out
}
