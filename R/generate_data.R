#' Sample from a spatial scale-GEV-model
#'
#' This function generates data with a max-stable-dependence structure and margins from
#' the scale-GEV-model.
#'
#' @param n Record length
#' @param d Number of stations (i.e. the dimension of the sampled data)
#' @param scale vector of length \eqn{d} containing the scale parameters for each station
#' @param loc vector of length \eqn{d} containing the location parameters for each station
#' @param shape vector of length \eqn{d} containing the shape parameters for each station
#' @param alpha vector of length \eqn{d} containing the trend parameters for each station
#' @param covmod The covariance model for the max-stable-process. Must be one of `gauss`,
#' `brown`,  `powexp`, `cauchy`, `whitmat`, `bessel`.
#' @param locations A \eqn{d \times 2} matrix containing longitude and latitude of the \eqn{d} stations.
#' @param temp.cov A vector of length \eqn{n} giving the temporal covariate.
#' @param MSparam A named list containing the parameters for the max-stable process. For `gauss`, these must
#' be `cov11, cov12, cov22`, for the other models these are `range, smooth`. For further details,
#' see \code{\link[SpatialExtremes]{rmaxstab}}.
#' @param seed This seed is set for reproducibility.
#'
#' @return A \eqn{n \times d} matrix containing containing the data.
#' @export
#'
#' @examples
#'
#' # with the same scale-GEV model parameters for d stations
#' plot.ts(generateData(seed = 1))
#' # with different scale-GEV model parameters at station 1 and 2
#' plot.ts(generateData(scale = c(2, 2, rep(5, 7)), shape =  c(0.2, 0.2, rep(0.1, 7)),  seed = 1))
generateData <- function(n = 100, d = 9,
                         scale = rep(5, d), loc = rep(20, d), shape = rep(0.1, d), alpha = rep(2.5, d),
                         covmod = "gauss", locations = matrix(20*abs(stats::rnorm(d*2)), ncol = 2),
                         temp.cov = (1:100)/100,
                         MSparam = list(cov11 = 0.4, cov12 = 0.2, cov22 = 0.9, range = 0.55, smooth = 1),
                         seed) {
  if(!(covmod %in% c("gauss","brown", "powexp", "cauchy", "whitmat", "bessel"))) {
    stop("Choose one of the following max-stable models: \n
         gauss, brown, powexp, cauchy, whitmat, bessel") }
  # generate data
  set.seed(seed)
  if(covmod == "gauss") {
    Xfr <- SpatialExtremes::rmaxstab(n = n,
                                     coord = locations,
                                     cov.mod = covmod, cov11 = MSparam$cov11,
                                     cov12 = MSparam$cov12,
                                     cov22 = MSparam$cov22)

  } else {

    Xfr <- SpatialExtremes::rmaxstab(n = n,
                                     coord = locations,
                                     cov.mod = covmod, nugget = 0,
                                     range = MSparam$range,
                                     smooth = MSparam$smooth)

  }
  Xgev <-  sapply(1:d, function(i){
    loc.i <- loc[i]*exp(alpha[i]*temp.cov/loc[i])
    scale.i <- scale[i]*exp(alpha[i]*temp.cov/loc[i])

    SpatialExtremes::frech2gev(Xfr[ , i], loc = loc.i, scale = scale.i,
                               shape = shape[i]) })



  colnames(Xgev) <- 1:d

  Xgev
}
