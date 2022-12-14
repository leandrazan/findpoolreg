% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/teststat_fun.R
\name{htheta}
\alias{htheta}
\title{Functionals to express Homogeneity/Equal Distribution hypotheses}
\usage{
htheta(par.mat, equal_distr = TRUE)
}
\arguments{
\item{par.mat}{A \eqn{4\times d} matrix containing the estimated parameter vectors for the \eqn{d} stations.}

\item{equal_distr}{Logical; if \code{TRUE} (default), the hypothesis is the one of equal distributions across the different
stations. If \code{FALSE}, the hypothesis is that the local scaling model holds across the region.}
}
\value{
A vector of length \eqn{4(d-1)} (equal distribution hypothesis) or \eqn{3(d-1)} (local scaling hypothesis).
}
\description{
Functionals \eqn{h} for which, under the null-hypothesis, we have \eqn{h(\theta)  = 0}.
}
\examples{
A <- matrix(1:8, ncol = 2)
htheta(A)
}
