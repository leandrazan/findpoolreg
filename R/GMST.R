#' Global Mean Surface Temperature anomaly w.r.t. 1951-1980
#'
#'
#' @format A data frame with columns
#' \describe{
#' \item{Year}{The year}
#' \item{annualGMST}{The global mean surface temperature anomaly of the
#' respective year, with respect to the mean of 1951-1980}
#' \item{smoothedGMST}{Global mean surface temperature anomaly smoothed over
#' 4 years. For the values of 2019-2021, these are smoothed version of 2016-2021
#'  values (because 2022+  values not observed the time of package construction.)}
#' }
#'
#' @source \url{https://data.giss.nasa.gov/gistemp/}
"GMST"
