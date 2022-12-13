# estimate 'regionwise return period/return level'


#' Generate one observation with a given dependence structure and given GEV
#' parameters
#'
#' @param ms_fit List containing the max-stable process parameters, i.e.
#' * `covmod`: A character string giving the covariance model of the max-stable process model
#'  (one of `gauss, brown, powexp, cauchy, whitmat, bessel`.)
#' * `parhat`: The parameters of the respective covariance model.
#' @param gev_par Values of the scale-GEV parameters, in the order \eqn{\mu, \sigma,\gamma, \alpha}.
#' @param coords The coordinates (longitude and latitude) of the locations for which you
#' want to generate data.
#' @param ref_gmst The reference value of the temporal covariate for which data is generated.
#' @param seed A random seed.
#'
#' @return A vector of length `nrow(coords)` containing one observation for each location.
#' @export
#'
#' @examples
#'
#' d <- 10
#' coords <- matrix(20*abs(stats::rnorm(d*2)), ncol = 2)
#' gen_gev_ms_sample(ms_fit = list(covmod = "gauss", parhat = c(1.2, 0.1, 1.2)),
#' gev_par = c(20, 4, 0.2, 2), coords = coords, ref_gmst = 0.9, seed = 1)
#'
gen_gev_ms_sample <- function(ms_fit, gev_par, coords, ref_gmst, seed) {
  set.seed(seed)

  d <- nrow(coords)
  if(is.list(ms_fit$parhat)) { ms_fit$parhat <- unlist(ms_fit$parhat)}
  if(ms_fit$covmod == "gauss") {
    Xfr <- SpatialExtremes::rmaxstab(n = 1, coord = coords, cov.mod = ms_fit$covmod,
                                     cov11 = ms_fit$parhat[1], cov12 = ms_fit$parhat[2],
                                     cov22 = ms_fit$parhat[3])
  }
  else if (ms_fit$covmod == "brown") {
    Xfr <- SpatialExtremes::rmaxstab(n = 1, coord = coords, cov.mod = ms_fit$covmod,
                                     range = ms_fit$parhat[1], smooth = ms_fit$parhat[2])
  }
  else {
    Xfr <- SpatialExtremes::rmaxstab(n = 1, coord = coords, cov.mod = ms_fit$covmod,
                                     nugget= ms_fit$parhat[1],  range = ms_fit$parhat[2], smooth = ms_fit$parhat[3])
  }

  mut <- gev_par[1]*exp(gev_par[4]*ref_gmst/gev_par[1])
  sigmat <- gev_par[2]*exp(gev_par[4]*ref_gmst/gev_par[1])

  Xgev <- sapply(1:d, function(i) {
    SpatialExtremes::frech2gev(Xfr[, i], loc = mut, scale = sigmat,
                               shape = gev_par[3])
  })

  names(Xgev) <- 1:d
  Xgev
}

#' Generate sample of independent regionwise maxima
#'
#' Calls the function `gen_gev_ms_sample` repeatedly and computes maximum of each
#' iteration.
#'
#' @inheritParams gen_gev_ms_sample
#' @param B The resulting sample length.
#' @param start Start value for the random seed; `gen_gev_ms_sample` is called with
#' seeds  `(start + 1) : (start + B)`.
#'
#' @return A vector of length `B` containing simulated regionwise maxima.
#' @export
#'
#' @examples
#'
#' d <- 10
#' coords <- matrix(20*abs(stats::rnorm(d*2)), ncol = 2)
#' maxis <- gen_sample_max(ms_fit = list(covmod = "gauss", parhat = c(1.2, 0.1, 1.2)),
#'    gev_par = c(20, 4, 0.2, 2), ref_gmst = 0.9, B = 100, coords = coords)
#' plot.ts(maxis)
gen_sample_max <- function(ms_fit, gev_par, ref_gmst, B, coords, start = 0) {


  purrr::map_dbl(1:B, ~ max(gen_gev_ms_sample(ms_fit = ms_fit, gev_par = gev_par, ref_gmst = ref_gmst,
                                              coords = coords, seed = start + .x)))

}

#' Evaluate parameters of scale-GEV at covariate values
#'
#' @param parvec Vector of scale-GEV parameters, in the order \eqn{\mu, \sigma,\gamma, \alpha}.
#' @param cvrt Vector of reference covariate values.
#'
#' @return  A vector of length 3 (when `cvrt` consists of only one reference value) or a
#' matrix of dimension 3 x `length(cvrt)`, containing values of
#' \eqn{\mu(t), \sigma(t), \gamma}.
#' @export
#'
#' @examples get_gev_par(c(20, 4, 00.2, 2), cvrt = c(0.9, 0.99))
get_gev_par <- function(parvec, cvrt) {

  n.cvrt <- length(cvrt)
  mu <- parvec[1]*exp(parvec[4]/parvec[1]*cvrt)
  sigma  <- parvec[2]*exp(parvec[4]/parvec[1]*cvrt)

  if(n.cvrt == 1) {
    return(c(mu, sigma, parvec[3]))
  }
  else {
    return(matrix(c(mu, sigma, rep(parvec[3], n.cvrt)), byrow = TRUE, nrow = 3))
  }

}


#' Estimate 'regionwise' RP and RL for pooled region
#'
#' @param data Data matrix containing data of the pooled region.
#' @param coords Coordinates of the locations that make up the pooled region.
#' @param temp.cov Numeric vector of same length `nrow(data)` containing the
#' temporal covariate.
#' @param ref_gmst Reference value of temporal covariate for which to compute the
#'  regionwise estimates.
#' @param B A large integer giving the sample size based on which the RP/RL is
#'  estimated.
#' @param r A numeric vector of values for which to compute the regionwise return
#'  period.
#' @param T.year  A numeric vèctor of values for which to compute the regionwise
#' return level.
#' @param ms_models A list of max-stable process models that are fitted to the
#' data in order to find the best fit.
#' @param return_ecdf Logical; whether to return the empirical cumulative distribution
#' function of the simulated samples.
#' @param start Passed to \code{\link[findpoolreg]{gen_sample_max}}, giving start
#' value for random seeds.
#'
#' @return A list containing a subset of the following components:
#' * `ecdf_max` The empirical cumulative distribution function of the simulated
#'  sample of regionwise maxima.
#'  * `rl_reg` The estimated regionwise RLs.
#'  * `rp_reg` The estimated regionwise RPs.
#'  * `gev_par` The estimated scale-GEV parameters of the pooled sample.
#'  * `fit_ms` A list containing the components `covmod` and `parhat`, giving
#'  the (covariance) model and the respective estimated parameters.
#'
#' @export
#'
#' @examples
#' d <- 10
#' coords <- matrix(20*abs(stats::rnorm(d*2)), ncol = 2)
#'
#'  xx <- generateData(locations = coords, d = 10, seed = 2)
#'
#' ## better use a larger value of B than done here in applications!
#' get_regionwise_estimates(data = xx, coords = coords, temp.cov = (1:100)/100,
#'     ref_gmst = 0.9, B = 5000, T.year = 100, r = 65)
get_regionwise_estimates <- function(data, coords, temp.cov, ref_gmst, B, r = NULL, T.year = NULL,
                                     ms_models = c("gauss", "brown", "whitmat", "powexp"),
                                     return_ecdf = TRUE, start = 0) {

  d.mat <- ncol(data)
  if(!(d.mat == nrow(coords))) { stop(" 'Data' must have as many columns as
                                           there are rows in 'coords'.")}

  if(is.null(r) & is.null(T.year) & !return_ecdf) { stop("Please provide either a value
                                                         for `r` or a value for `T.year` (or both).")}

  gev_fit_pool <- findpoolreg::fit_scalegev(data, temp.cov = temp.cov)

  for (i in 1:d.mat) {
    data[, i] <- findpoolreg::scalegev2frech(data = data[, i], temp.cov = temp.cov,
                                                par = gev_fit_pool$mle)
  }

  fitted_ms <- list()
  for (ms.mod in ms_models) {
    fit.tmp <- tryCatch(SpatialExtremes::fitmaxstab(data,
                                                    coord = coords, ms.mod, warn = FALSE), error = function(a) {
                                                      list(fitted.values = NA)
                                                    })
    fitted_ms <- if (!anyNA(fit.tmp$fitted.values)) {
      append(fitted_ms, list(fit.tmp))
    }
  }

  ms_model <- fitted_ms[[which.min(purrr::map_dbl(fitted_ms, ~ SpatialExtremes::TIC(.x)))]]
  ms_fit_sub <- list()
  ms_fit_sub$parhat <- ms_model$fitted.values
  ms_fit_sub$covmod <- ms_model$cov.mod
  if(ms_fit_sub$covmod == "Gaussian"){ ms_fit_sub$covmod <- "gauss" }

  maxis <- findpoolreg::gen_sample_max(ms_fit = ms_fit_sub,
                          gev_par = gev_fit_pool$mle, ref_gmst = ref_gmst, B = B, coords = coords, start = start)

  ecdf_max <- stats::ecdf(maxis)

  if(return_ecdf){
    list_ecdf <- list("ecdf_max" = ecdf_max)
  }
  else {
    list_ecdf <- list()
  }

  if(!is.null(T.year)) {
    rl_reg <- stats::quantile(maxis, p = 1-1/T.year)
    names(rl_reg) <- paste("T =", T.year)
    list_ecdf <- append(list_ecdf, list(rl_reg = rl_reg))
  }
  if(!is.null(r)) {
    rp_reg <-  1/(1-ecdf_max(r))
    names(rp_reg) <- paste( "r =", r)
    list_ecdf <- append(list_ecdf,  list(rp_reg = rp_reg))
  }

  append(list_ecdf,  list(gev_par = gev_fit_pool$mle,
                          fit_ms = ms_fit_sub) )


}


