% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loglik_fun.R
\name{nll_scalegev}
\alias{nll_scalegev}
\title{Negative log-likelihood of univariate GEV-scale-model}
\usage{
nll_scalegev(params, data, temp.cov)
}
\arguments{
\item{params}{A vector of length 4}

\item{data}{The vector of observations}

\item{temp.cov}{Vector containing values of the temporal covariate.}
}
\value{
Numeric value of the log-likelihood of the data evaluated at the given parameters.
}
\description{
Negative log-likelihood of univariate GEV-scale-model
}
\examples{
# generate some data
xx <- exp((1:100)/100)*evd::rgev(100)

nll_scalegev(params = c(2,1, 0, 1), data = xx, temp.cov = (1:100)/100)
}
