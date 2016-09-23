% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bSpline.R
\name{bSpline}
\alias{bSpline}
\title{B-Spline Basis for Polynomial Splines}
\usage{
bSpline(x, df = NULL, knots = NULL, degree = 3, intercept = FALSE,
        Boundary.knots = range(x), ...)
}
\arguments{
\item{x}{The predictor variable.  Missing values are allowed and will be
returned as they were.}

\item{df}{Degrees of freedom.  One can specify \code{df} rather than
\code{knots}, then the function chooses "df - degree"
(minus one if there is an intercept) knots at suitable quantiles of \code{x}
(which will ignore missing values).  The default, \code{NULL}, corresponds
to no inner knots, i.e., "degree - intercept".}

\item{knots}{The internal breakpoints that define the spline.  The
default is \code{NULL}, which results in a basis for ordinary
polynomial regression.  Typical values are the mean or median
for one knot, quantiles for more knots.  See also
\code{Boundary.knots}.}

\item{degree}{Non-negative integer degree of the piecewise polynomial. The
default value is 3 for cubic splines. Zero degree is allowed for this
function, which is the only difference compared with
\code{\link[splines]{bs}} in package \code{splines}.}

\item{intercept}{If \code{TRUE}, an intercept is included in the basis;
Default is \code{FALSE}.}

\item{Boundary.knots}{Boundary points at which to anchor the B-spline basis.
By default, they are the range of the non-\code{NA} data.  If both
\code{knots} and \code{Boundary.knots} are supplied, the basis parameters
do not depend on \code{x}. Data can extend beyond \code{Boundary.knots}.}

\item{...}{Optional arguments for future usage.}
}
\value{
A matrix of dimension \code{length(x)} by
\code{df = degree + length(knots)} (plus one if intercept is included).
Attributes that correspond to the arguments specified are returned
for usage for \code{\link{predict.bSpline2}}.
}
\description{
This function generates the B-spline basis matrix for a polynomial spline.
}
\details{
It is an augmented function of \code{\link[splines]{bs}} in package
\code{splines} for B-spline basis that allows piecewise constant spline
basis with zero degree. When the argument \code{degree} is greater than
zero, it internally calls \code{\link[splines]{bs}} and generates a basis
matrix for representing the family of piecewise polynomials with the
specified interior knots and degree, evaluated at the values of \code{x}.
The function has the same arguments with \code{\link[splines]{bs}} for ease
usage.
}
\examples{
library(graphics)
x <- seq(0, 1, by = 0.01)
knots <- c(0.3, 0.5, 0.6)
bMat <- bSpline(x, knots = knots, degree = 0, intercept = TRUE)
matplot(x, bMat, type = "l", ylab = "B-spline basis with degree zero")
abline(v = knots, lty = 2, col = "gray")
}
\seealso{
\code{\link{predict.bSpline2}} for evaluation at given (new) values;
\code{\link{ibs}} for integral of B-spline basis;
\code{\link{mSpline}} for M-spline basis;
\code{\link{iSpline}} for I-spline basis;
\code{\link{cSpline}} for C-spline basis.
}
