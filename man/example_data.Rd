% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/example_data.R
\docType{data}
\name{example_data}
\alias{example_data}
\title{Example data set}
\format{
A 75 x 9 matrix of observations sampled from a scale-GEV model with spatial dependence structure
described by a max-stable process.
}
\usage{
example_data
}
\description{
Example data set
}
\examples{
# The data was generated with the following code.
meancoords <-  as.matrix(findpoolreg::example_grid[ , 6:7])
cvrt <- findpoolreg::GMST$smoothedGMST[68:142]
generateData(n = 75, d = 9,
            scale = c(6, 6, 6, 5, 5, 6, 5, 5, 6),
            loc =  c(23, 23, 22, 20, 20, 20, 20, 20, 21),
            shape =  c( 0.12, 0.1, 0.15, 0.1, 0.1, 0.15, 0.11, 0.1, 0.05),
            alpha =  c(rep(3, 3), rep(2.5, 6)),
            seed = 1, locations = meancoords, temp.cov = cvrt)

}
\keyword{datasets}
