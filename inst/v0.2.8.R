##' B-Spline Basis for Polynomial Splines
##'
##' This function generates the B-spline basis matrix for a polynomial spline.
##'
##' It is an augmented function of \code{\link[splines]{bs}} in package
##' \code{splines} for B-spline basis that allows piecewise constant (close on
##' the left, open on the right) spline basis with zero degree. When the
##' argument \code{degree} is greater than zero, it internally calls
##' \code{\link[splines]{bs}} and generates a basis matrix for representing the
##' family of piecewise polynomials with the specified interior knots and
##' degree, evaluated at the values of \code{x}.  The function has the same
##' arguments with \code{\link[splines]{bs}} for ease usage.
##'
##' @usage
##' bSpline(x, df = NULL, knots = NULL, degree = 3L, intercept = FALSE,
##'         Boundary.knots = range(x, na.rm = TRUE), ...)
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
##' @examples
##' library(splines2)
##' x <- seq.int(0, 1, 0.01)
##' knots <- c(0.3, 0.5, 0.6)
##' bsMat <- bSpline(x, knots = knots, degree = 0, intercept = TRUE)
##'
##' library(graphics)
##' matplot(x, bsMat, type = "l", ylab = "Piecewise constant B-spline bases")
##' abline(v = knots, lty = 2, col = "gray")
##' @seealso
##' \code{\link{predict.bSpline2}} for evaluation at given (new) values;
##' \code{\link{dbs}}, \code{\link{deriv.bSpline2}} for derivatives;
##' \code{\link{ibs}} for integral of B-splines;
##' \code{\link{mSpline}} for M-splines;
##' \code{\link{iSpline}} for I-splines;
##' \code{\link{cSpline}} for C-splines.
##' @importFrom splines bs
##' @importFrom stats stepfun
##' @export
bSpline <- function(x, df = NULL, knots = NULL, degree = 3L, intercept = FALSE,
                    Boundary.knots = range(x, na.rm = TRUE), ...)
{
    ## check and reformat 'degree'
    if ((degree <- as.integer(degree)) < 0)
        stop("'degree' must be a nonnegative integer.")

    ## sort and remove possible NA's in internal knots if exist
    if (length(knots))
        knots <- sort.int(knots)

    ## take care of possible NA's in `x`
    nax <- is.na(x)
    if (all(nax))
        stop("The 'x' cannot be all NA's!")

    ## call splines::bs for non-zero degree
    if (degree > 0) {
        out <- splines::bs(x = x, df = df, knots = knots,
                           degree = degree, intercept = intercept,
                           Boundary.knots = Boundary.knots)
        ## add "x" to attributes
        attr(out, "x") <- x
        ## throw out warning if any internal knot outside boundary.knots
        knots <- attr(out, "knots")
        Boundary.knots <- attr(out, "Boundary.knots")
        ## any internal knots placed outside of boundary knots?
        outside_knots <- (knots <= Boundary.knots[1L]) |
            (knots >= Boundary.knots[2L])
        if (any(outside_knots))
            warning(wrapMessages(
                "Some internal knots were not placed",
                "inside of boundary knots,",
                "which may cause \nill-conditioned bases!"
            ))
        ## update classes
        class(out) <- c("matrix", "bSpline2")
        return(out)
    }

    ## else degree is zero
    ## remove NA's in x
    xx <- if (nas <- any(nax)) x[! nax] else x

    ## check whether any of x is outside of the boundary knots
    outside_x <- rep(FALSE, length(xx))
    if (! missing(Boundary.knots)) {
        if (! is.numeric(Boundary.knots) || anyNA(Boundary.knots))
            stop(wrapMessages(
                "The 'Boundary.knots' has to be",
                "numeric vector of length 2",
                "with no missing value."
            ))
        if (length(Boundary.knots) > 2) {
            warning(wrapMessages(
                "Only the first two values",
                "in the 'Boundary.knots' were used."
            ))
            Boundary.knots <- Boundary.knots[seq_len(2L)]
        }
        Boundary.knots <- sort.int(Boundary.knots)
        outside_x <- (xx < Boundary.knots[1L]) | (xx > Boundary.knots[2L])
    }
    if (any(outside_x))
        warning(wrapMessages(
            "Some 'x' values beyond boundary knots",
            "may cause ill-conditioned bases!"
        ))

    ## prepare inputs for piecewise constant bases
    inputs <- pieceConst(x = xx[! outside_x],
                         df = df, knots = knots,
                         Boundary.knots = Boundary.knots,
                         intercept = intercept)
    knots <- inputs$knots
    ## potentially, df is a bad name since df is also a function in stats
    df <- inputs$df

    ## piecewise constant basis
    augKnots <- c(Boundary.knots[1L], knots,
                  Boundary.knots[2L] + 10 * .Machine$double.eps)
    bsMat <- sapply(seq_len(df), function (i) {
        foo <- stats::stepfun(augKnots[i: (i + 1L)], c(0L, 1L, 0L))
        foo(xx)
    })

    ## close on the right boundary knot for the last constant piece?
    ## if (any(rightX <- xx == Boundary.knots[2L]))
    ##     bsMat[rightX, df] <- 1

    ## make sure bsMat is a matrix
    if (! is.matrix(bsMat))
        bsMat <- matrix(bsMat, nrow = length(xx))

    ## include intercept or not
    if (! intercept) {
        bsMat <- bsMat[, - 1L, drop = FALSE]
    }

    ## keep NA's as is
    if (nas) {
        nmat <- matrix(NA, length(nax), ncol(bsMat))
        nmat[! nax, ] <- bsMat
        bsMat <- nmat
    }

    ## add dimnames for consistency with bs returns
    row.names(bsMat) <- names(x)
    colnames(bsMat) <- as.character(seq_len(df - as.integer(! intercept)))

    ## on attributes
    tmp <- list(degree = degree,
                knots = if (is.null(knots)) numeric(0L) else knots,
                Boundary.knots = Boundary.knots,
                intercept = intercept, x = x)
    attributes(bsMat) <- c(attributes(bsMat), tmp)
    class(bsMat) <- c("matrix", "bSpline2")
    bsMat
}


### internal function ==========================================================
##' @importFrom stats quantile
pieceConst <- function (x, df, knots, Boundary.knots, intercept)
{
    ind <- (is.null(df) + 1L) * is.null(knots) + 1L
    ## ind == 1: knots is not NULL;
    ##   df0 = df = length(knots) + 1L

    ## ind == 2: df is not NULL, while knots is NULL;
    ## df := function input
    ##   number of knots = df - as.integer(intercept) from `splines::bs`
    ## df0 := DF of spline bases from splines definition
    ##      = length(knots) + 1L
    ##      = df - as.integer(intercept) + 1L

    ## ind == 3: both df and knots are NULL; one-piece constant
    ##   number of knots = 0, df0 = 1

    int_intercept <- as.integer(intercept)
    df0 <- switch(ind,
                  length(knots) + 1L,
                  {
                      int_df <- as.integer(df)
                      if (int_df < 1L)
                          stop("The spepcified `df` must be positive!",
                               call. = FALSE)
                      int_df - int_intercept + 1L
                  },
                  {
                      if (! intercept)
                          stop(wrapMessages(
                              "The 'intercept' has to be 'TRUE'",
                              "for one-piece const basis."
                          ), call. = FALSE)
                      1L
                  })
    if (ind > 1L) {
        tknots <- df0 + 1L
        quans <- seq.int(from = 0, to = 1,
                         length.out = tknots)[- c(1L, tknots)]
        knots <- as.numeric(stats::quantile(x, quans))
    } else {
        ## any internal knots placed outside of boundary knots?
        outside_knots <- (knots <= Boundary.knots[1L]) |
            (knots >= Boundary.knots[2L])
        ## remove internal knots placed outside of boundary knots
        if (any(outside_knots)) {
            knots <- knots[! outside_knots]
            df0 <- df0 - sum(outside_knots)
            warning(wrapMessages(
                "Only internal knots placed inside of",
                "the boundary knots were considered."
            ), call. = FALSE)
        }
        if (! is.null(df) && df != df0)
            warning(wrapMessages(
                "The 'df' specified was not appropriate.",
                sprintf("Used 'df = %d' instead.", df0)
            ),  call. = FALSE)
    }
    list(df = df0, knots = knots)
}


##' C-Spline Basis for Polynomial Splines
##'
##' This function generates the convex regression spline (called C-spline) basis
##' matrix by integrating I-spline basis for a polynomial spline.
##'
##' It is an implementation of the close form C-spline basis derived from
##' the recursion formula of I-spline and M-spline.  Internally, it calls
##' \code{\link{iSpline}} and generates a basis matrix for representing the
##' family of piecewise polynomials and their corresponding integrals with the
##' specified interior knots and degree, evaluated at the values of \code{x}.
##'
##' @usage
##' cSpline(x, df = NULL, knots = NULL, degree = 3L, intercept = TRUE,
##'         Boundary.knots = range(x, na.rm = TRUE), scale = TRUE, ...)
##'
##' @param x The predictor variable.  Missing values are allowed and will be
##'     returned as they were.
##' @param df Degrees of freedom.  One can specify \code{df} rather than
##'     \code{knots}, then the function chooses "df - degree" (minus one if
##'     there is an intercept) knots at suitable quantiles of \code{x} (which
##'     will ignore missing values).  The default, \code{NULL}, corresponds to
##'     no inner knots, i.e., "degree - intercept".
##' @param knots The internal breakpoints that define the spline.  The default
##'     is \code{NULL}, which results in a basis for ordinary polynomial
##'     regression.  Typical values are the mean or median for one knot,
##'     quantiles for more knots.  See also \code{Boundary.knots}.
##' @param degree Non-negative integer degree of the piecewise polynomial. The
##'     default value is 3 for cubic splines.
##' @param intercept If \code{TRUE} by default, all spline bases are included.
##'     Notice that when using C-Spline for shape-restricted regression,
##'     \code{intercept = TRUE} should be set even when an intercept term is
##'     considered additional to the spline bases in the model.
##' @param Boundary.knots Boundary points at which to anchor the C-spline basis.
##'     By default, they are the range of the non-\code{NA} data.  If both
##'     \code{knots} and \code{Boundary.knots} are supplied, the basis
##'     parameters do not depend on \code{x}. Data can extend beyond
##'     \code{Boundary.knots}.
##' @param scale Logical value (\code{TRUE} by default) indicating whether
##'     scaling on C-spline basis is required. If TRUE, C-spline basis is scaled
##'     to have unit height at right boundary knot; the corresponding I-spline
##'     and M-spline basis matrices shipped in attributes are also scaled to the
##'     same extent.
##' @param ... Optional arguments for future usage.
##'
##' @return A matrix of dimension \code{length(x)} by
##' \code{df = degree + length(knots)} (plus on if intercept is included).
##' The attributes that correspond to the arguments specified are returned
##' for the usage of other functions in this package.
##' @references
##' Meyer, M. C. (2008). Inference using shape-restricted regression splines.
##' \emph{The Annals of Applied Statistics}, 1013--1033. Chicago
##' @examples
##' library(splines2)
##' x <- seq.int(0, 1, 0.01)
##' knots <- c(0.3, 0.5, 0.6)
##'
##' ### when 'scale = TRUE' (by default)
##' csMat <- cSpline(x, knots = knots, degree = 2)
##'
##' library(graphics)
##' matplot(x, csMat, type = "l", ylab = "C-spline basis")
##' abline(v = knots, lty = 2, col = "gray")
##' isMat <- deriv(csMat)
##' msMat <- deriv(csMat, derivs = 2)
##' matplot(x, isMat, type = "l", ylab = "scaled I-spline basis")
##' matplot(x, msMat, type = "l", ylab = "scaled M-spline basis")
##'
##' ### when 'scale = FALSE'
##' csMat <- cSpline(x, knots = knots, degree = 2, scale = FALSE)
##' ## the corresponding I-splines and M-splines (with same arguments)
##' isMat <- iSpline(x, knots = knots, degree = 2)
##' msMat <- mSpline(x, knots = knots, degree = 2, intercept = TRUE)
##' ## or using deriv methods (more efficient)
##' isMat1 <- deriv(csMat)
##' msMat1 <- deriv(csMat, derivs = 2)
##' ## equivalent
##' stopifnot(all.equal(isMat, isMat1, check.attributes = FALSE))
##' stopifnot(all.equal(msMat, msMat1, check.attributes = FALSE))
##' @seealso
##' \code{\link{predict.cSpline}} for evaluation at given (new) values;
##' \code{\link{deriv.cSpline}} for derivatives;
##' \code{\link{iSpline}} for I-splines;
##' \code{\link{mSpline}} for M-splines.
##' @importFrom stats stepfun
##' @export
cSpline <- function(x, df = NULL, knots = NULL, degree = 3L, intercept = TRUE,
                    Boundary.knots = range(x, na.rm = TRUE), scale = TRUE, ...)
{
    ## I-spline basis for inputs
    isOut <- iSpline(x = x, df = df, knots = knots, degree = degree,
                     intercept = intercept, Boundary.knots = Boundary.knots)

    ## update input
    degree <- attr(isOut, "degree")
    knots <- attr(isOut, "knots")
    bKnots <- attr(isOut, "Boundary.knots")
    ord <- 1L + degree
    nKnots <- length(knots)
    df <- nKnots + ord

    ## take care of possible NA's in `x` for the following calculation
    nax <- is.na(x)
    if ((nas <- any(nax)))
        x <- x[! nax]
    nX <- length(x)

    ## define knot sequence
    aKnots <- sort(c(rep(bKnots, ord + 1L), knots))

    ## generate I-spline basis with (degree + 1)
    augX <- c(x, bKnots[2L])
    isOut1 <- iSpline(x = augX, knots = knots, degree = ord,
                      intercept = FALSE, Boundary.knots = bKnots)

    ## function determining j from x
    j <- if (length(knots)) {
             foo <- stats::stepfun(x = knots, y = seq.int(ord, df))
             as.integer(foo(augX))
         } else {
             rep.int(ord, nX + 1L)
         }

    numer1 <- diff(aKnots, lag = ord + 1L)[- 1L]
    ## if there is at least one internal knot
    if (nKnots) {
        ## calculate C-spline basis at each internal knot t_j
        isOutKnots <- iSpline(knots, knots = knots, degree = ord,
                              intercept = FALSE, Boundary.knots = bKnots)
        matKnots <- rep(numer1, each = nKnots) * isOutKnots / (ord + 1)
        augKnots <- seq_len(nKnots) + ord
        diffKnots <- diff(knots)
        csKnots <- lapply(seq_len(nKnots), function(i, idx) {
            ji <- augKnots[i]
            a <- matKnots[i, ]
            js <- seq_len(ji)
            a[- js] <- 0
            a[js] <- rev(cumsum(rev(a[js])))
            a[idx < ji - ord] <- diffKnots[ji - ord - 1L]
            a
        }, idx = seq_len(df))
        csKnots <- do.call(rbind, csKnots)

        idxMat <- lower.tri(csKnots, diag = TRUE)
        linList <- lapply(seq_len(nKnots), function(ind) {
            cumsum(csKnots[idxMat[, ind], ind])
        })
        csKnots[idxMat] <- do.call(c, linList)
    } else {
        csKnots <- matrix(0, 1L, df)
    }

    ## calculate C-spline basis at each x
    matX <- rep(numer1, each = nX + 1) * isOut1 / (ord + 1)
    csOut <- lapply(seq_len(nX + 1L), function(i, idx) {
        ji <- j[i]
        xx <- augX[i]
        a <- matX[i, ]
        js <- seq_len(ji)
        a[- js] <- 0
        a[js] <- rev(cumsum(rev(a[js])))
        a[idx < ji - ord] <- xx - knots[ji - ord] +
            csKnots[ji - ord, idx < ji - ord]
        a
    }, idx = seq_len(df))
    csOut <- do.call(rbind, csOut)

    if (! intercept)
        csOut <- csOut[, - 1L, drop = FALSE]
    scl <- unname(csOut[nX + 1L, ])
    csOut <- csOut[- (nX + 1L), ]

    ## mSpline basis matrix
    msMat <- attr(isOut, "msMat")

    ## keep NA's as is for csOut
    if (nas) {
        nmat <- matrix(NA, length(nax), ncol(csOut))
        nmat[! nax, ] <- csOut
        csOut <- nmat
    }

    ## scale C-spline, I-spline, and M-spline basis
    if (scale) {
        vec <- rep(1 / scl, each = length(nax))
        csOut <- vec * csOut
        isOut <- vec * isOut
        msMat <- vec * msMat
        attr(isOut, "scale") <- attr(msMat, "scale") <- scale
        attr(isOut, "scales") <- attr(msMat, "scales") <- scl
    }

    ## output
    attr(isOut, "msMat") <- NULL
    attributes(csOut) <- c(attributes(isOut),
                           list(isMat = isOut, msMat = msMat,
                                scale = scale, scales = scl))
    attr(csOut, "derivs") <- NULL
    class(csOut) <- c("matrix", "cSpline")
    csOut
}


##' Derivative of B-Spline Basis for Polynomial Splines
##'
##' This function produces the derivative of given order of B-splines.  It is an
##' implementation of the close form derivative of B-spline basis based on
##' recursion relation.  At knots, the derivative is defined to be the right
##' derivative.
##'
##' The function is similar with \code{\link[splines]{splineDesign}}. However,
##' it provides a more user-friendly interface, a more considerate \code{NA}'s
##' handling.  Internally, it calls \code{\link{bSpline}} and generates a basis
##' matrix for representing the family of piecewise polynomials and their
##' corresponding derivative with the specified interior knots and degree,
##' evaluated at the values of \code{x}. The function \code{splineDesign} in
##' \code{splines} package can also be used to calculate derivative of
##' B-splines.
##'
##' @usage
##' dbs(x, derivs = 1L, df = NULL, knots = NULL, degree = 3L,
##'     intercept = FALSE, Boundary.knots = range(x, na.rm = TRUE), ...)
##'
##' @param x The predictor variable.  Missing values are allowed and will be
##'     kept and returned as they were.
##' @param derivs A positive integer specifying the order of derivative.  By
##'     default, it is \code{1L} for the first derivative.
##' @param df Degrees of freedom of the B-spline basis to be differentiated.
##'     One can specify \code{df} rather than \code{knots}, then the function
##'     chooses "df - degree" (minus one if there is an intercept) knots at
##'     suitable quantiles of \code{x} (which will ignore missing values).  The
##'     default, \code{NULL}, corresponds to no inner knots, i.e.,
##'     "degree - intercept".
##' @param knots The internal breakpoints that define the B-spline basis to be
##'     differentiated.  The default is \code{NULL}, which results in a basis
##'     for ordinary polynomial regression.  Typical values are the mean or
##'     median for one knot, quantiles for more knots.  See also
##'     \code{Boundary.knots}.
##' @param degree Non-negative integer degree of the piecewise polynomial to be
##'     differentiated. The default value is 3 for the integral of cubic
##'     B-splines.
##' @param intercept If \code{TRUE}, an intercept is included in the basis;
##'     Default is \code{FALSE}.
##' @param Boundary.knots Boundary points at which to anchor the B-spline basis
##'     to be differentiated. By default, they are the range of the
##'     non-\code{NA} data.  If both \code{knots} and \code{Boundary.knots} are
##'     supplied, the basis parameters do not depend on \code{x}.
##' @param ... Optional arguments for future usage.
##'
##' @return A matrix of dimension \code{length(x)} by
##' \code{df = degree + length(knots)} (plus on if intercept is included).
##' Attributes that correspond to the arguments specified are returned
##' for usage of other functions in this package.
##' @references
##' De Boor, Carl. (1978). \emph{A practical guide to splines}.
##' Vol. 27. New York: Springer-Verlag.
##' @examples
##' library(splines2)
##' x <- seq.int(0, 1, 0.01)
##' knots <- c(0.2, 0.4, 0.7)
##' ## the second derivative of cubic B-splines with three internal knots
##' dMat <- dbs(x, derivs = 2L, knots = knots, intercept = TRUE)
##'
##' ## compare with the results from splineDesign
##' ord <- attr(dMat, "degree") + 1L
##' bKnots <- attr(dMat, "Boundary.knots")
##' aKnots <- c(rep(bKnots[1L], ord), knots, rep(bKnots[2L], ord))
##' res <- splines::splineDesign(aKnots, x = x, derivs = 2L)
##' stopifnot(all.equal(res, dMat, check.attributes = FALSE))
##' @seealso
##' \code{\link{predict.dbs}} for evaluation at given (new) values;
##' \code{\link{deriv.dbs}} for derivative method;
##' \code{\link{bSpline}} for B-splines;
##' \code{\link{ibs}} for integral of B-splines.
##' @importFrom stats quantile
##' @export
dbs <- function(x, derivs = 1L, df = NULL, knots = NULL, degree = 3L,
                intercept = FALSE, Boundary.knots = range(x, na.rm = TRUE), ...)
{
    ## check order of derivative
    derivs <- as.integer(derivs)
    if (derivs < 1L)
        stop("'derivs' has to be a positive integer.")

    ## check and reformat 'degree'
    if ((degree <- as.integer(degree)) < 0)
        stop("'degree' must be a nonnegative integer.")

    ## sort and remove possible NA's in internal knots if exist
    if (length(knots))
        knots <- sort.int(knots)

    ## take care of possible NA's in `x`
    nax <- is.na(x)
    if (all(nax))
        stop("'x' cannot be all NA's!")
    nas <- any(nax)
    ## remove NA's
    xx <- if (nas) x[! nax] else x

    ## check Boundary.knots specified by users
    outside <- rep(FALSE, length(xx))
    if (! missing(Boundary.knots)) {
        Boundary.knots <- sort(Boundary.knots[seq_len(2)])
        outside <- (xx < Boundary.knots[1L]) | (xx > Boundary.knots[2L])
    }

    ## determine knots from df if missing
    inter <- as.integer(intercept)
    if (! is.null(df)) {
        df0 <- length(knots) + degree + inter
        if (tmp <- (df < df0))
            warning(sprintf("'df' was too small; have used %d", df0))

        df <- ifelse(tmp, df0, df)
        nKnots <- df - degree - inter
        if (is.null(knots) && nKnots > 0) {
            quans <- seq.int(from = 0, to = 1,
                             length.out = nKnots + 2L)[- c(1L, nKnots + 2L)]
            knots <- stats::quantile(xx[! outside], quans)
        }
    }
    ## update degree of freedom from inputs
    df0 <- length(knots) + degree + 1L
    df <- df0 - 1L + inter

    ## attribute knots for output
    knotsAttr <- if (is.null(knots)) numeric(0L) else knots

    ## for derivs > degree
    if (derivs > degree) {
        ## df == 0, i.e., no basis returned
        if (! df)
            warning("Degree of freedom is zero.")
        dMat <- matrix(0, nrow = length(x), ncol = df)
        if (nas)
            dMat[nax, ] <- NA
        tmp <- list(degree = degree,
                    knots = knotsAttr,
                    Boundary.knots = Boundary.knots,
                    intercept = intercept,
                    x = x, derivs = derivs)
        attributes(dMat) <- c(attributes(dMat), tmp)
        class(dMat) <- c("matrix", "dbs")
        return(dMat)
    }

    ## B-spline bases
    dMat <- bSpline(xx, knots = knots, degree = degree - derivs,
                    intercept = TRUE, Boundary.knots = Boundary.knots, ...)

    ## derivative matrix
    for (iter in seq_len(derivs)) {
        ## define knot sequence according to the bases being differentiated
        ord <- degree - derivs + iter + 1L
        aKnots <- sort(c(rep(Boundary.knots, ord), knots))
        denom <- diff(aKnots, lag = ord - 1L)
        facVec <- ifelse(denom > 0, (ord - 1L) / denom, 0)
        dMat0 <- cbind(0, dMat, 0)
        dMat <- sapply(seq_len(df0 - derivs + iter), function(a)
        {
            idx <- a : (a + 1L)
            tmpMat <- dMat0[, idx, drop = FALSE]
            facVec[idx[1L]] * tmpMat[, 1L, drop = FALSE] -
                facVec[idx[2L]] * tmpMat[, 2L, drop = FALSE]
        })
        ## recover dimension after sapply
        if (! is.matrix(dMat))
            dMat <- matrix(dMat, nrow = 1L)
    }

    ## take care of intercept
    if (! intercept)
        dMat <- dMat[, - 1L, drop = FALSE]

    ## keep NA's as is
    if (nas) {
        nmat <- matrix(NA, length(nax), df)
        nmat[! nax, ] <- dMat
        dMat <- nmat
    }

    ## add dimnames for consistency with returns from splines::bs
    row.names(dMat) <- names(x)
    colnames(dMat) <- as.character(seq_len(df))

    ## on attributes
    tmp <- list(degree = degree,
                knots = knotsAttr,
                Boundary.knots = Boundary.knots,
                intercept = intercept,
                x = x, derivs = derivs)
    attributes(dMat) <- c(attributes(dMat), tmp)
    class(dMat) <- c("matrix", "dbs")
    ## return
    dMat
}


##' Derivative of Splines
##'
##' \code{deriv} methods that obtain derivative of given order of B-splines,
##' M-spline, I-splines, and C-splines, etc. At knots, the derivative is defined
##' to be the right derivative. By default, the function returns the first
##' derivative. For derivatives of order greater than one, the nested call such
##' as \code{deriv(deriv(expr))} is supported but not recommended. For a better
##' performance, argument \code{derivs} should be specified instead.
##'
##' The function is designed for most of the objects generated from this
##' package. It internally extracts necessary information about the input spline
##' basis matrix from its attributes. So the function will not work if some
##' attribute is not available.
##'
##' @name deriv
##'
##' @param expr Objects of class \code{bSpline2}, \code{ibs}, \code{dbs},
##'     \code{mSpline}, \code{iSpline}, or \code{cSpline}, etc.
##' @param derivs A positive integer specifying the order of derivatives. By
##'     default, it is \code{1L} for the first derivative.
##' @param ... Other arguments for further usage.
##'
##' @return A matrix of dimension \code{length(x)} by
##' \code{df = degree + length(knots)} (plus one if intercept is included).
##' Attributes that correspond to the arguments specified are returned
##' for usage for other function in this package.
##' @references
##' De Boor, Carl. (1978). \emph{A practical guide to splines}.
##' Vol. 27. New York: Springer-Verlag.
##' @examples
##' library(splines2)
##' x <- c(seq.int(0, 1, 0.1), NA)  # NA's will be kept.
##' knots <- c(0.3, 0.5, 0.6)
##'
##' ## integal of B-splines and the corresponding B-splines integrated
##' ibsMat <- ibs(x, knots = knots)
##' bsMat <- bSpline(x, knots = knots)
##'
##' ## the first derivative
##' d1Mat <- deriv(ibsMat)
##' stopifnot(all.equal(bsMat, d1Mat, check.attributes = FALSE))
##'
##' ## the second derivative
##' d2Mat1 <- deriv(bsMat)
##' d2Mat2 <- deriv(ibsMat, derivs = 2L)
##' ## nested calls are supported but not recommended
##' d2Mat3 <- deriv(deriv(ibsMat))
##' stopifnot(all.equal(d2Mat1, d2Mat2, d2Mat3, check.attributes = FALSE))
##'
##' ## C-splines, I-splines, M-splines and the derivatives
##' csMat <- cSpline(x, knots = knots, scale = FALSE)
##' isMat <- iSpline(x, knots = knots)
##' stopifnot(all.equal(isMat, deriv(csMat), check.attributes = FALSE))
##'
##' msMat <- mSpline(x, knots = knots, intercept = TRUE)
##' stopifnot(all.equal(msMat, deriv(isMat), deriv(csMat, 2),
##'                     deriv(deriv(csMat)), check.attributes = FALSE))
##'
##' dmsMat <- mSpline(x, knots = knots, intercept = TRUE, derivs = 1)
##' stopifnot(all.equal(dmsMat, deriv(msMat), deriv(isMat, 2),
##'                     deriv(deriv(isMat)), deriv(csMat, 3),
##'                     deriv(deriv(deriv(csMat))), check.attributes = FALSE))
##' @seealso
##' \code{\link{bSpline}} for B-splines;
##' \code{\link{ibs}} for integral of B-splines;
##' \code{\link{mSpline}} for M-splines;
##' \code{\link{iSpline}} for I-splines;
##' \code{\link{cSpline}} for C-splines.
##' @importFrom stats deriv
NULL


##' @rdname deriv
##' @export
deriv.bSpline2 <- function(expr, derivs = 1L, ...)
{
    attr(expr, "derivs") <- derivs
    dMat <- do.call(dbs, attributes(expr))
    class(dMat) <- c("matrix", "dbs")
    dMat
}


##' @rdname deriv
##' @export
deriv.dbs <- function(expr, derivs = 1L, ...)
{
    attr(expr, "derivs") <- attr(expr, "derivs") + derivs
    dMat <- do.call(dbs, attributes(expr))
    class(dMat) <- c("matrix", "dbs")
    dMat
}


##' @rdname deriv
##' @export
deriv.ibs <- function(expr, derivs = 1L, ...)
{
    ## quick check on derivs
    derivs <- as.integer(derivs)
    if (derivs < 1L)
        stop("'derivs' has to be a positive integer.")

    ## if first derivative, take result from existing attribute
    if (derivs == 1L) {
        out <- attr(expr, "bsMat")
        if (is.null(out))
            out <- do.call(bSpline, attributes(expr))
        return(out)
    }

    ## for higher order of derivative
    attr(expr, "derivs") <- derivs - 1L
    dMat <- do.call(dbs, attributes(expr))
    class(dMat) <- c("matrix", "dbs")
    dMat
}


##' @rdname deriv
##' @export
deriv.mSpline <- function(expr, derivs = 1L, ...)
{
    ## call function dbs
    derivs0 <- attr(expr, "derivs")
    attr(expr, "derivs") <- ifelse(is.null(derivs0), derivs, derivs0 + derivs)
    dMat <- do.call(mSpline, attributes(expr))

    ## for possible scaling of objects from deriv.cSpline
    scale <- attr(expr, "scale")
    scl <- attr(expr, "scales")
    if (isTRUE(scale)) {
        dMat <- dMat / rep(scl, each = nrow(dMat))
        attr(dMat, "scale") <- scale
        attr(dMat, "scales") <- scl
    }
    ## prepare for output
    class(dMat) <- c("matrix", "mSpline")
    dMat
}



##' @rdname deriv
##' @export
deriv.iSpline <- function(expr, derivs = 1L, ...)
{
    ## quick check on derivs
    derivs <- as.integer(derivs)
    if (derivs < 1L)
        stop("'derivs' has to be a positive integer.")

    ## extract existing result from attributes for the first derivative
    if (derivs == 1L) {
        dMat <- attr(expr, "msMat")
        if (is.null(dMat))
            dMat <- deriv.mSpline(expr, derivs = 1L, ...)
    } else {
        ## for derivative of higher order
        dMat <- deriv.mSpline(expr = expr, derivs = derivs - 1L, ...)
    }

    class(dMat) <- c("matrix", "mSpline")
    dMat
}


##' @rdname deriv
##' @export
deriv.cSpline <- function(expr, derivs = 1L, ...)
{
    derivs <- as.integer(derivs)
    if (derivs < 1L)
        stop("'derivs' has to be a positive integer.")

    scale <- attr(expr, "scale")
    scl <- attr(expr, "scales")
    if (derivs == 1L) {
        dMat <- attr(expr, "isMat")
        attr(dMat, "msMat") <- attr(expr, "msMat")
        class(dMat) <- c("matrix", "iSpline")
    } else if (derivs == 2L) {
        dMat <- attr(expr, "msMat")
        class(dMat) <- c("matrix", "mSpline")
    } else {
        dMat <- deriv.mSpline(expr = expr, derivs = derivs - 2L, ...)
        attr(dMat, "derivs") <- derivs - 2L
    }

    dMat
}

##' Integral of B-Spline Basis for Polynomial Splines
##'
##' This function generates the integral of B-spline basis matrix
##' for a polynomial spline. The arguments are exactly the same with function
##' \code{\link[splines]{bs}} in package \code{splines}.
##'
##' It is an implementation of the close form integral of B-spline basis based
##' on recursion relation.  Internally, it calls \code{\link{bSpline}} and
##' generates a basis matrix for representing the family of piecewise
##' polynomials and their corresponding integrals with the specified interior
##' knots and degree, evaluated at the values of \code{x}.
##'
##' @usage
##' ibs(x, df = NULL, knots = NULL, degree = 3, intercept = FALSE,
##'     Boundary.knots = range(x, na.rm = TRUE), ...)
##'
##' @param x The predictor variable.  Missing values are allowed and will be
##'     returned as they were.
##' @param df Degrees of freedom of the B-spline basis to be integrated.  One
##'     can specify \code{df} rather than \code{knots}, then the function
##'     chooses "df - degree" (minus one if there is an intercept) knots at
##'     suitable quantiles of \code{x} (which will ignore missing values).  The
##'     default, \code{NULL}, corresponds to no inner knots, i.e.,
##'     "degree - intercept".
##' @param knots The internal breakpoints that define the B-spline basis to be
##'     integrated.  The default is \code{NULL}, which results in a basis for
##'     ordinary polynomial regression.  Typical values are the mean or median
##'     for one knot, quantiles for more knots.  See also \code{Boundary.knots}.
##' @param degree Non-negative integer degree of the piecewise polynomial to be
##'     integrated. The default value is 3 for the integral of cubic B-splines.
##' @param intercept If \code{TRUE}, an intercept is included in the basis;
##'     Default is \code{FALSE}.
##' @param Boundary.knots Boundary points at which to anchor the B-spline basis
##'     to be integrated. By default, they are the range of the non-\code{NA}
##'     data.  If both \code{knots} and \code{Boundary.knots} are supplied, the
##'     basis parameters do not depend on \code{x}. Data can extend beyond
##'     \code{Boundary.knots}.
##' @param ... Optional arguments for future usage.
##'
##' @return A matrix of dimension \code{length(x)} by
##' \code{df = degree + length(knots)} (plus on if intercept is included).
##' Attributes that correspond to the arguments specified are returned
##' for usage of other functions in this package.
##' @references
##' De Boor, Carl. (1978). \emph{A practical guide to splines}.
##' Vol. 27. New York: Springer-Verlag.
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
##' @seealso
##' \code{\link{predict.ibs}} for evaluation at given (new) values;
##' \code{\link{deriv.ibs}} for derivative method.
##' \code{\link{bSpline}} for B-splines;
##' \code{\link{dbs}} for derivatives of B-splines;
##' @export
ibs <- function(x, df = NULL, knots = NULL, degree = 3, intercept = FALSE,
                Boundary.knots = range(x, na.rm = TRUE), ...)
{
    ## B-spline basis for inputs
    bsOut <- bSpline(x = x, df = df, knots = knots, degree = degree,
                     intercept = intercept, Boundary.knots = Boundary.knots)

    ## update input
    degree <- attr(bsOut, "degree")
    knots <- attr(bsOut, "knots")
    bKnots <- attr(bsOut, "Boundary.knots")
    ord <- 1L + degree

    ## define knot sequence
    aKnots <- sort(c(rep(bKnots, ord), knots))

    ## generate B-spline basis with (degree + 1)
    bsOut1 <- bSpline(x = x, knots = knots, degree = ord,
                      intercept = FALSE, Boundary.knots = bKnots)
    numer1 <- diff(aKnots, lag = ord)
    if (! intercept) {
        bsOut1 <- bsOut1[, - 1L, drop = FALSE]
        numer1 <- numer1[- 1L]
    }
    numer2 <- apply(bsOut1, 1, function(a) rev(cumsum(rev(a))))
    ibsOut <- t(numer1 * numer2) / ord

    ## output
    attributes(ibsOut) <- c(attributes(bsOut),
                            list(bsMat = bsOut, x = x))
    class(ibsOut) <- c("matrix", "ibs")
    ibsOut
}


##' I-Spline Basis for Polynomial Splines or its derivatives
##'
##' This function generates the I-spline (integral of M-spline) basis matrix for
##' a polynomial spline or its derivatives of given order..
##'
##' It is an implementation of the close form I-spline basis based on the
##' recursion formula of B-spline basis.  Internally, it calls
##' \code{\link{mSpline}} and \code{\link{bSpline}}, and generates a basis
##' matrix for representing the family of piecewise polynomials and their
##' corresponding integrals with the specified interior knots and degree,
##' evaluated at the values of \code{x}.
##'
##' @usage
##' iSpline(x, df = NULL, knots = NULL, degree = 3L, intercept = TRUE,
##'         Boundary.knots = range(x, na.rm = TRUE), derivs = 0L, ...)
##'
##' @param x The predictor variable.  Missing values are allowed and will be
##'     returned as they were.
##' @param df Degrees of freedom.  One can specify \code{df} rather than
##'     \code{knots}, then the function chooses "df - degree" (minus one if
##'     there is an intercept) knots at suitable quantiles of \code{x} (which
##'     will ignore missing values).  The default, \code{NULL}, corresponds to
##'     no inner knots, i.e., "degree - intercept".
##' @param knots The internal breakpoints that define the spline.  The default
##'     is \code{NULL}, which results in a basis for ordinary polynomial
##'     regression.  Typical values are the mean or median for one knot,
##'     quantiles for more knots.  See also \code{Boundary.knots}.
##' @param degree Non-negative integer degree of the piecewise polynomial. The
##'     default value is 3 for cubic splines. Note that the degree of I-spline
##'     is defined to be the degree of the associated M-spline instead of actual
##'     polynomial degree. In other words, I-spline basis of degree 2 is defined
##'     as the integral of associated M-spline basis of degree 2.
##' @param intercept If \code{TRUE} by default, all spline bases are included.
##'     Notice that when using I-Spline for monotonic regression,
##'     \code{intercept = TRUE} should be set even when an intercept term is
##'     considered additional to the spline bases in the model.
##' @param Boundary.knots Boundary points at which to anchor the I-spline basis.
##'     By default, they are the range of the non-\code{NA} data.  If both
##'     \code{knots} and \code{Boundary.knots} are supplied, the basis
##'     parameters do not depend on \code{x}. Data can extend beyond
##'     \code{Boundary.knots}.
##' @param derivs A non-negative integer specifying the order of derivatives of
##'     I-splines.
##' @param ... Optional arguments for future usage.
##'
##' @return A matrix of dimension \code{length(x)} by
##' \code{df = degree + length(knots)} (plus on if intercept is included).
##' Attributes that correspond to the arguments specified are returned
##' for usage of other functions in this package.
##' @references
##' Ramsay, J. O. (1988). Monotone regression splines in action.
##' \emph{Statistical science}, 3(4), 425--441.
##' @examples
##' ## Example given in the reference paper by Ramsay (1988)
##' library(splines2)
##' x <- seq.int(0, 1, by = 0.01)
##' knots <- c(0.3, 0.5, 0.6)
##' isMat <- iSpline(x, knots = knots, degree = 2)
##'
##' library(graphics)
##' matplot(x, isMat, type = "l", ylab = "I-spline basis")
##' abline(v = knots, lty = 2, col = "gray")
##'
##' ## the derivative of I-splines is M-spline
##' msMat1 <- iSpline(x, knots = knots, degree = 2, derivs = 1)
##' msMat2 <- mSpline(x, knots = knots, degree = 2, intercept = TRUE)
##' stopifnot(all.equal(msMat1, msMat2))
##' @seealso
##' \code{\link{predict.iSpline}} for evaluation at given (new) values;
##' \code{\link{deriv.iSpline}} for derivative method;
##' \code{\link{mSpline}} for M-splines;
##' \code{\link{cSpline}} for C-splines;
##' @importFrom stats stepfun
##' @export
iSpline <- function(x, df = NULL, knots = NULL, degree = 3L, intercept = TRUE,
                    Boundary.knots = range(x, na.rm = TRUE), derivs = 0L, ...)
{
    ## check order of derivative
    if (! missing(derivs)) {
        derivs <- as.integer(derivs)
        if (derivs < 0L)
            stop("'derivs' has to be a non-negative integer.")
    }

    ## M-spline basis for outputs in attributes
    msOut <- mSpline(x = x, df = df, knots = knots,
                     degree = degree, intercept = intercept,
                     Boundary.knots = Boundary.knots, derivs = 0L, ...)

    ## update input
    degree <- attr(msOut, "degree")
    knots <- attr(msOut, "knots")
    bKnots <- attr(msOut, "Boundary.knots")
    ord <- 1L + degree
    nKnots <- length(knots)
    df <- nKnots + ord

    ## default, for derivs == 0L, return I-splines
    if (! derivs) {
        ## define knot sequence
        ## aKnots <- sort(c(rep(bKnots, ord + 1L), knots))

        ## take care of possible NA's in `x` for the following calculation
        nax <- is.na(x)
        if (nas <- any(nax))
            x <- x[! nax]

        ## function determining j from x
        j <- if (nKnots) {
                 foo <- stats::stepfun(x = knots, y = seq.int(ord, df))
                 as.integer(foo(x))
             } else {
                 rep.int(ord, length(x))
             }

        ## calculate I-spline basis at non-NA x's
        ## directly based on B-spline
        bsOut1 <- bSpline(x = x, knots = knots, degree = ord,
                          intercept = FALSE, Boundary.knots = bKnots)

        isOut <- lapply(seq_along(j), function(i, idx) {
            a <- bsOut1[i, ]
            js <- seq_len(j[i])
            a[- js] <- 0
            a[js] <- rev(cumsum(rev(a[js])))
            a[idx < j[i] - ord] <- 1        # <=> a[idx < j[i] - degree] <- 1
            a
        }, idx = seq_len(df))
        isOut <- do.call(rbind, isOut)

        ## Or based on M-spline
        ## generate M-spline basis with (degree + 1)

        ## msOut1 <- mSpline(x = x, knots = knots, degree = ord,
        ##                   intercept = FALSE, Boundary.knots = bKnots)
        ## df <- length(knots) + ord
        ## numer1 <- diff(aKnots, lag = ord + 1)[- 1L]
        ## msMat <- rep(numer1, each = length(x)) * msOut1 / (ord + 1)
        ## msAugMat <- cbind(j, msMat)
        ## isOut <- t(apply(msAugMat, 1, function(b, idx = seq_len(df)) {
        ##     j <- b[1L]
        ##     a <- b[- 1L]
        ##     js <- seq_len(j)
        ##     a[- js] <- 0
        ##     a[js] <- rev(cumsum(rev(a[js])))
        ##     a[idx < j - ord] <- 1            # <=> a[idx < j - degree] <- 1
        ##     a
        ## }))

        ## intercept
        if (! intercept)
            isOut <- isOut[, - 1L, drop = FALSE]

        ## keep NA's as is
        if (nas) {
            nmat <- matrix(NA, length(nax), ncol(isOut))
            nmat[! nax, ] <- isOut
            isOut <- nmat
        }

    } else {
        ## for derivatives >= 1L
        out <- mSpline(x = x, df = df, knots = knots,
                       degree = degree, intercept = intercept,
                       Boundary.knots = Boundary.knots,
                       derivs = derivs - 1L, ...)
        return(out)
    }
    ## output
    attributes(isOut) <- c(attributes(msOut), list(msMat = msOut))
    class(isOut) <- c("matrix", "iSpline")
    isOut
}


### some trivial internal functions ============================================
## wrap messages and keep proper line length
wrapMessages <- function(..., strwrap.args = list()) {
    x <- paste(...)
    wrap_x <- do.call(strwrap, c(list(x = x), strwrap.args))
    paste(wrap_x, collapse = "\n")
}

## warning if x contains NA (or NaN)
na_warning <- function(x, sub_env = c("current", "parent", "grandparent"),
                       num_grandparent = 2L, ...)
{
    sub_env <- switch(
        match.arg(sub_env),
        "current" = environment(),
        "parent" = parent.frame(),
        "grandparent" = parent.frame(num_grandparent)
    )
    objName = deparse(substitute(x, sub_env))
    if (anyNA(x))
        warning(wrapMessages(
            sprintf("Found `NA` values in `%s`.", objName)
        ), call. = FALSE)
    invisible(x)
}

## stop if x contains NA (or NaN)
na_stop <- function(x, sub_env = c("current", "parent", "grandparent"),
                       num_grandparent = 2L, ...)
{
    sub_env <- switch(
        match.arg(sub_env),
        "current" = environment(),
        "parent" = parent.frame(),
        "grandparent" = parent.frame(num_grandparent)
    )
    objName = deparse(substitute(x, sub_env))
    if (anyNA(x))
        stop(wrapMessages(
            sprintf("Found `NA` values in `%s`.", objName)
        ), call. = FALSE)
    invisible(x)
}

## is x a numeric matrix (optionally of nRow rows and nCol columns)
isNumMatrix <- function(x, nRow = NULL, nCol = NULL,
                        warn_na = TRUE, error_na = ! warn_na,
                        sub_env = "parent", ...)
{
    out <- is.numeric(x) && is.matrix(x)
    if (out) {
        nDim <- dim(x)
        if (! is.null(nRow)) out <- out && nDim[1L] == nRow
        if (! is.null(nCol)) out <- out && nDim[2L] == nCol
        if (error_na) na_stop(x, sub_env = sub_env, ...)
        if (warn_na) na_warning(x, sub_env = sub_env, ...)
    }
    out
}


##' M-Spline Basis for Polynomial Splines and its Derivatives
##'
##' This function generates the basis matrix of the regression spline called
##' M-spline or its derivatives of given order.  For monotone regression,
##' \code{\link{iSpline}} should be used.
##'
##' It is an implementation of the close form M-spline basis based on
##' relationship between M-spline basis and B-spline basis.  In fact, M-spline
##' basis is a rescaled version of B-spline basis. Internally, it calls function
##' \code{\link{bSpline}} and generates a basis matrix for representing the
##' family of piecewise polynomials with the specified interior knots and
##' degree, evaluated at the values of \code{x}.
##'
##' @usage
##' mSpline(x, df = NULL, knots = NULL, degree = 3L, intercept = FALSE,
##'         Boundary.knots = range(x, na.rm = TRUE), derivs = 0L, ...)
##'
##' @param x The predictor variable.  Missing values are allowed and will be
##'     returned as they were.
##' @param df Degrees of freedom.  One can specify \code{df} rather than
##'     \code{knots}, then the function chooses "df - degree" (minus one if
##'     there is an intercept) knots at suitable quantiles of \code{x} (which
##'     will ignore missing values).  The default, \code{NULL}, corresponds to
##'     no inner knots, i.e., "degree - intercept".
##' @param knots The internal breakpoints that define the spline.  The default
##'     is \code{NULL}, which results in a basis for ordinary polynomial
##'     regression.  Typical values are the mean or median for one knot,
##'     quantiles for more knots.  See also \code{Boundary.knots}.
##' @param degree Non-negative integer degree of the piecewise polynomial. The
##'     default value is 3 for cubic splines. Zero degree is allowed for
##'     piecewise constant basis.
##' @param intercept If \code{TRUE}, all bases will be returned. The default
##'     value is \code{FALSE}.
##' @param Boundary.knots Boundary points at which to anchor the M-spline basis.
##'     By default, they are the range of the non-\code{NA} data.  If both
##'     \code{knots} and \code{Boundary.knots} are supplied, the basis
##'     parameters do not depend on \code{x}. Data can extend beyond
##'     \code{Boundary.knots}.
##' @param derivs A non-negative integer specifying the order of derivatives of
##'     M-splines. The default value is \code{0L} for M-spline bases.
##' @param ... Optional arguments for future usage.
##'
##' @return A matrix of dimension \code{length(x)} by
##' \code{df = degree + length(knots)} (plus one if intercept is included).
##' Attributes that correspond to the arguments specified are returned
##' for usage of other functions in this package.
##' @references
##' Ramsay, J. O. (1988). Monotone regression splines in action.
##' \emph{Statistical science}, 3(4), 425--441.
##' @examples
##' ## Example given in the reference paper by Ramsay (1988)
##' library(splines2)
##' x <- seq.int(0, 1, 0.01)
##' knots <- c(0.3, 0.5, 0.6)
##' msMat <- mSpline(x, knots = knots, degree = 2, intercept = TRUE)
##'
##' library(graphics)
##' matplot(x, msMat, type = "l", ylab = "M-spline basis")
##' abline(v = knots, lty = 2, col = "gray")
##'
##' ## derivatives of M-splines
##' dmsMat <- mSpline(x, knots = knots, degree = 2,
##'                   intercept = TRUE, derivs = 1)
##' ## or using the 'deriv' method
##' dmsMat1 <- deriv(msMat)
##' stopifnot(all.equal(dmsMat, dmsMat1, check.attributes = FALSE))
##' @seealso
##' \code{\link{predict.mSpline}} for evaluation at given (new) values;
##' \code{\link{deriv.mSpline}} for derivative method;
##' \code{\link{bSpline}} for B-splines;
##' \code{\link{iSpline}} for I-splines;
##' \code{\link{cSpline}} for C-splines.
##' @export
mSpline <- function(x, df = NULL, knots = NULL, degree = 3L, intercept = FALSE,
                    Boundary.knots = range(x, na.rm = TRUE), derivs = 0L, ...)
{
    ## check order of derivative
    if (! missing(derivs)) {
        derivs <- as.integer(derivs)
        if (derivs < 0L)
            stop("'derivs' has to be a non-negative integer.")
    }

    bsOut <- if (derivs) {
                 dbs(x = x, derivs = derivs, df = df, knots = knots,
                     degree = degree, intercept = intercept,
                     Boundary.knots = Boundary.knots, ...)
             } else {
                 bSpline(x = x, df = df, knots = knots,
                         degree = degree, intercept = intercept,
                         Boundary.knots = Boundary.knots, ...)
             }

    ## update input
    ord <- attr(bsOut, "degree") + 1L
    knots <- attr(bsOut, "knots")
    bKnots <- attr(bsOut, "Boundary.knots")

    ## define knot sequence
    aKnots <- sort(c(rep(bKnots, ord), knots))

    ## transformation from B-splines to M-splines
    denom <- diff(aKnots, lag = ord)
    transCoef <- ifelse(denom > 0, ord / denom, 0)
    if (! intercept)
        transCoef <- transCoef[- 1L]
    msOut <- rep(transCoef, each = length(x)) * bsOut
    attr(msOut, "derivs") <- derivs
    class(msOut) <- c("matrix", "mSpline")
    msOut
}


##' Evaluate a Spline Basis
##'
##' This function evaluates a predefined spline basis at (new) given values.
##'
##' These are methods for the generic function \code{predict} for objects
##' inheriting from class \code{bSpline2}, \code{ibs}, \code{mSpline},
##' \code{iSpline}, or \code{cSpline}.  If \code{newx} is not given, the
##' function returns the input object.  For object returned by function
##' \code{\link{cSpline}}, the \code{mSpline} and \code{iSpline} objects shipped
##' in attributes should not be evaluated by this function if \code{rescale} is
##' \code{TRUE}.  See \code{\link{cSpline}} for details.
##'
##' @name predict
##' @param object Objects of class \code{bSpline2}, \code{ibs}, \code{mSpline},
##'     \code{iSpline}, or \code{cSpline} having attributes describing
##'     \code{knots}, \code{degree}, etc.
##' @param newx The \code{x} values at which evaluations are required.
##' @param ... Optional argument for future usage.
##'
##' @return An object just like the \code{object} input, except evaluated at
##' the new values of \code{x}.
##'
##' @examples
##' library(splines2)
##' x <- seq.int(0, 1, 0.2)
##' knots <- c(0.3, 0.5, 0.6)
##' newX <- seq.int(0.1, 0.9, 0.2)
##'
##' ## for B-splines
##' bsMat <- bSpline(x, knots = knots, degree = 2)
##' predict(bsMat, newX)
##'
##' ## for integral of B-splines
##' ibsMat <- ibs(x, knots = knots, degree = 2)
##' predict(ibsMat, newX)
##'
##' ## for derivative of B-splines
##' dbsMat <- dbs(x, knots = knots, degree = 2)
##' predict(dbsMat, newX)
##'
##' ## for M-spline
##' msMat <- mSpline(x, knots = knots, degree = 2)
##' predict(msMat, newX)
##'
##' ## for I-spline
##' isMat <- iSpline(x, knots = knots, degree = 2)
##' predict(isMat, newX)
##'
##' ## for C-spline
##' csMat <- cSpline(x, knots = knots, degree = 2)
##' predict(csMat, newX)
##' @seealso
##' \code{\link{bSpline}} for B-splines;
##' \code{\link{ibs}} for integral of B-splines;
##' \code{\link{dbs}} for derivative of B-splines;
##' \code{\link{mSpline}} for M-splines;
##' \code{\link{iSpline}} for I-splines;
##' \code{\link{cSpline}} for C-splines.
##' @importFrom stats predict
NULL


##' @rdname predict
##' @export
predict.bSpline2 <- function(object, newx, ...)
{
    if (missing(newx))
        return(object)
    a <- c(list(x = newx),
           attributes(object)[c("degree", "knots", "Boundary.knots",
                                "intercept")])
    do.call("bSpline", a)
}


##' @rdname predict
##' @export
predict.ibs <- function(object, newx, ...)
{
    if (missing(newx))
        return(object)
    a <- c(list(x = newx),
           attributes(object)[c("degree", "knots", "Boundary.knots",
                                "intercept")])
    do.call("ibs", a)
}


##' @rdname predict
##' @export
predict.dbs <- function(object, newx, ...)
{
    if (missing(newx))
        return(object)
    a <- c(list(x = newx),
           attributes(object)[c("degree", "knots", "Boundary.knots",
                                "intercept", "derivs")])
    do.call("dbs", a)
}


##' @rdname predict
##' @export
predict.mSpline <- function(object, newx, ...)
{
    if (missing(newx))
        return(object)
    a <- c(list(x = newx),
           attributes(object)[c("degree", "knots", "Boundary.knots",
                                "intercept", "derivs")])
    do.call("mSpline", a)
}


##' @rdname predict
##' @export
predict.iSpline <- function(object, newx, ...)
{
    if (missing(newx))
        return(object)
    a <- c(list(x = newx),
           attributes(object)[c("degree", "knots", "Boundary.knots",
                                "intercept", "derivs")])
    do.call("iSpline", a)
}


##' @rdname predict
##' @export
predict.cSpline <- function(object, newx, ...)
{
    if (missing(newx))
        return(object)
    a <- c(list(x = newx),
           attributes(object)[c("degree", "knots", "Boundary.knots",
                                "intercept", "scale")])
    do.call("cSpline", a)
}

##' Print Out a Spline Basis Matrix
##'
##' \code{Print} methods that simply print out the spline basis matrix without
##' unnecessary attributes.
##'
##' @name print
##' @param x Objects of class \code{bSpline2}, \code{ibs}, \code{mSpline},
##' \code{iSpline}, or \code{cSpline}, etc.
##' @param ... Optional argument for future usage.
##'
##' @return Object input.
NULL


##' @rdname print
##' @export
print.bSpline2 <- function(x, ...) {
    print.default(tidyAttr(x, ...))
    invisible(x)
}


##' @rdname print
##' @export
print.ibs <- function(x, ...) {
    print.default(tidyAttr(x, ...))
    invisible(x)
}


##' @rdname print
##' @export
print.dbs <- function(x, ...) {
    print.default(tidyAttr(x, ...))
    invisible(x)
}


##' @rdname print
##' @export
print.mSpline <- function(x, ...) {
    print.default(tidyAttr(x, ...))
    invisible(x)
}

##' @rdname print
##' @export
print.iSpline <- function(x, ...) {
    print.default(tidyAttr(x, ...))
    invisible(x)
}


##' @rdname print
##' @export
print.cSpline <- function(x, ...) {
    print.default(tidyAttr(x, ...))
    invisible(x)
}

## remove all attributes but dim and dimnames
tidyAttr <- function(x, ...) {
    dimen <- attr(x, "dim")
    dimenName <- attr(x, "dimnames")
    attributes(x) <- NULL
    attr(x, "dim") <- dimen
    attr(x, "dimnames") <- dimenName
    x
}
