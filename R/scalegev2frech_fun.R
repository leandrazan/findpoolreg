#' Transform scale-GEV margins to unit Frechet margins
#'
#' @param data Numeric vector of observations
#' @param temp.cov Numeric vector of temporal covariate, with same length as `data`
#' @param par Either `NULL` (the default) or a vector containing the location, scale, shape and
#' trend parameter of the scale-GEV model. If `NULL`, the parameters are estimated with
#' maximum likelihood.
#'
#' @return A vector with same length as input data.
#' @export
#'
#' @examples
#' xx <- generateData(n = 100, temp.cov = (1:100)/100, seed=1)
#'
#' # apply transformation on first column
#' xx[, 1] <- scalegev2frech(xx[, 1], temp.cov = (1:100)/100)
#'
#' # spot the difference
#' plot.ts(xx)
scalegev2frech <- function(data, temp.cov, par = NULL) {


  if(!is.null(ncol(data))) { stop("Data must be univariate. Use spat_scalegev2frech() instead.")}
  if(!(length(data) == length(temp.cov))) { stop("Data and temporal covariate must have same length.")}

  if(is.null(par)) {
    par <- fit_scalegev(data, temp.cov = temp.cov, hessian = FALSE)$mle
  }

  loc <- par[1]*exp(par[4]*temp.cov/par[1])
  scale <- par[2]*exp(par[4]*temp.cov/par[1])
  unitfrech <- gev2frech(data, loc = loc, scale = scale,
                                          shape = par[3])
  unitfrech
}

#' Transform scale-GEV margins to unit Frechet margins for multivariate data
#' @param data Numeric matrix of observations. Each column corresponds to one station.
#' @param temp.cov Numeric vector of temporal covariate, with same length as `nrow(data)`
#' @param par Either `NULL` (the default) or a \eqn{4 \times ncol(data)} matrix containing the location, scale, shape and
#' trend parameters of the scale-GEV model of each column of the data. If `NULL`, the parameters are estimated stationwise with
#' maximum likelihood.
#'
#' @return A matrix of same dimension as input data.
#' @export
#'
#' @examples
#' xx <- generateData(n = 100, temp.cov = (1:100)/100, seed=1)
#'
#' # apply transformation
#' yy <- spat_scalegev2frech(xx, temp.cov = (1:100)/100)
#'
#' # spot the difference
#' plot.ts(xx)
#' plot.ts(yy)
spat_scalegev2frech <- function(data, temp.cov, par = NULL) {

  d <- ncol(data)
  if(is.null(d)) { stop("Dimension of data must be at least 2. Use scalegev2frech() instead.")}
  for( i in 1:d) {
    data[, i] <- scalegev2frech(data = data[, i], temp.cov = temp.cov, par = par[, i])
  }
  data
}
