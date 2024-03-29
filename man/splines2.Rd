% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splines2-package.R
\docType{package}
\name{splines2}
\alias{splines2}
\alias{splines2-package}
\title{splines2: Regression Spline Functions and Classes}
\description{
This package provides functions to construct basis matrices of
\itemize{
\item B-splines
\item M-splines
\item I-splines
\item convex splines (C-splines)
\item periodic splines
\item natural cubic splines
\item generalized Bernstein polynomials
\item along with their integrals (except C-splines) and derivatives
    of given order by closed-form recursive formulas
}
}
\details{
In addition to the R interface, it also provides a C++ header-only library
integrated with \pkg{Rcpp}, which allows the construction of spline basis
functions directly in C++ with the help of \pkg{Rcpp} and
\pkg{RcppArmadillo}.  Thus, it can also be treated as one of the \pkg{Rcpp*}
packages.  A toy example package that uses the C++ interface is available at
<https://github.com/wenjie2wang/example-pkg-Rcpp-splines2>.

The package \pkg{splines2} is intended to be a user-friendly supplement to
the base package \pkg{splines}.  The trailing number two in the package name
means "too" (and by no means refers to the generation two).  See Wang and
Yan (2021) for details and illustrations of how the package can be applied
to shape-restricted regression.
}
\references{
Wang, W., & Yan, J. (2021). Shape-restricted regression splines with R
package \pkg{splines2}. \emph{Journal of Data Science}, 19(3), 498--517.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://wwenjie.org/splines2}
  \item \url{https://github.com/wenjie2wang/splines2}
  \item Report bugs at \url{https://github.com/wenjie2wang/splines2/issues}
}

}
\author{
\strong{Maintainer}: Wenjie Wang \email{wang@wwenjie.org} (\href{https://orcid.org/0000-0003-0363-3180}{ORCID})

Authors:
\itemize{
  \item Jun Yan (\href{https://orcid.org/0000-0003-4401-7296}{ORCID})
}

}
\keyword{internal}
