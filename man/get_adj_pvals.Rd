% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrap_fun.R
\name{get_adj_pvals}
\alias{get_adj_pvals}
\title{Compute adjusted p-values}
\usage{
get_adj_pvals(
  tibres,
  methods = c("holm", "BY", "BH"),
  rejection = FALSE,
  level = 0.1
)
}
\arguments{
\item{tibres}{A tibble containing a column named \code{p_boot}, as returned by \code{bootstrap_scalegev_subsets}.
Can also be a vector of p-values.}

\item{methods}{A subset of \verb{holm, BY, BH}, giving the method that is used to adjust the p-values.}

\item{rejection}{Logical; whether to add an additional column giving the rejection result at the given level.}

\item{level}{Rejection level; only needed when \code{rejection = TRUE}.}
}
\value{
The original tibble extended with the adjusted p-values (and rejection results when \code{rejection = TRUE}).
Columns containing the adjusted p-values are named with \verb{p_holm, p_bh, p_by}.
}
\description{
Compute adjusted p-values
}
\examples{

\dontrun{
 # generate 6-dimensional data where two stations deviate from the rest
cvrt <- (1:100)/100
coords <- matrix(20*abs(stats::rnorm(6*2)), ncol = 2)
x <- generateData(seed = 2, n = 100, temp.cov = cvrt, d= 6, locations = coords,
loc = c(rep(10,4), c(8, 8)), scale = c(rep(2,4), c(1.5, 1.5)), shape = rep(0.1, 6),
alpha = rep(2,6))

sbsts <- list( c(1,2), c(1,3), c(1,4), c(1,5), c(1,6))
bootres <- bootstrap_subsets_ms(data = x, temp.cov = cvrt,
             locations = coords, varmeth = "chain", B = 200, subsets = sbsts)
get_adj_pvals(bootres) # H0 ist rejected at 5\%-level for the deviating stations 5 and 6
}
}
