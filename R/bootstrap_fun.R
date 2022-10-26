

#' Generate bootstrap samples with dependence structure from data and unit frechet margins
#'
#' This function transform margins to unit Frechet, fits several max-stable models to
#' the transformed data and chooses the best fit. Then ootstrap samples from this fitted max-stable model
#' are generated.
#'
#' @param data Numeric matrix of observations. Each column corresponds to one station.
#' @param temp.cov  Numeric vector of temporal covariate, with same length as `nrow(data)`
#' @param locations A matrix containing longitude and latitude of the stations. Each row corresponds to one station.
#' @param ms_models A vector containing the names of the max-stable models to fit. Must be a subset of
#' `gauss, brown, powexp, cauchy, whitmat, bessel`.
#' @param B The number of bootstrap samples that are generated.
#' @param sel_crit The selection criterion based on which the max-stable model is selected. Must be one of
#' `TIC, AIC`.
#' @param warn_msfit logical; whether to print warnings possibly generated during fitting the max-stable process.
#'
#' @return Returns a list of  `B` bootstrap samples, each being a matrix of same dimension as input data.
#' @export
#'
#' @examples
#' # generate data
#'
#' d <- 8
#' coords <-  matrix(20*abs(stats::rnorm(d*2)), ncol = 2)
#' xx <- generateData(seed = 1, d = d, locations = coords,
#'  scale = c(1:d), loc = rep(20, d), shape = seq(0, 0.3, length.out = d), alpha = rep(2.5, d))
#'
#' # generate bootstrap samples
#' bootsamps <- generate_bootsamp_unitfrech(xx, temp.cov = (1:100)/100, locations = coords,
#'                                         ms_models = c("powexp", "gauss", "brown"), B = 2,
#'                                          sel_crit = "TIC", warn_msfit = FALSE)
#' plot.ts(bootsamps[[1]])
#' plot.ts(bootsamps[[2]])
generate_bootsamp_unitfrech <- function(data, temp.cov, locations,
                                        ms_models = c("powexp", "gauss", "brown"), B = 500, sel_crit = "TIC",
                                        warn_msfit = TRUE) {

  if(!all(ms_models %in% c("gauss","brown", "powexp", "cauchy", "whitmat", "bessel"))) {
    stop("Only implemented for max-stable models in \n
    'gauss', 'brown', 'powexp', 'cauchy', 'whitmat', 'bessel'.")
  }

  if(!(sel_crit %in% c("TIC", "AIC"))) { stop("Selection criterion must be one of 'TIC', 'AIC'.")}

  ## transform margins to unit Frechet based on local GEV parameters
  unitfrech <- spat_scalegev2frech(data = data, temp.cov = temp.cov)

  ## Fit max-stable processes using models specified in ms_models (margins are unit Frechet)
  fitted_ms <- list()
  for( ms.mod in ms_models) {
    fitted_ms <- append(fitted_ms, list(SpatialExtremes::fitmaxstab(unitfrech, coord = locations, ms.mod, warn = warn_msfit)))
  }

  ## choose best fit based on selection criterion
  if(sel_crit == "TIC") {
    sel.vals <- purrr::map_dbl(fitted_ms, ~ SpatialExtremes::TIC(.x))
    if(all(is.na(sel.vals))) {
      sel_crit <- "AIC"
      warning("Selection criterion was changed to AIC because TIC didn't exist.")}
  }

  if(sel_crit == "AIC") {
    sel.vals <- purrr::map_dbl(fitted_ms, ~ AIC(.x))
  }

  model.sel <- which.min(sel.vals)

  # chosen fitted max-stable model
  fitted <- fitted_ms[[model.sel]]

  # covariance model of chosen max-stable process
  Kovmod <- fitted$cov.mod
  if(Kovmod == "Gaussian") { Kovmod <- "gauss"}

  ## estimated parameters of the max-stable process
  parhat <- fitted$fitted.values

  # generate bootstrap samples
  if(Kovmod == "gauss"){
    cov11_hat <- parhat[1]
    cov12_hat <- parhat[2]
    cov22_hat <- parhat[3]

    # sim data
    # NAMatrix <- is.na(data)

    X_star <- replicate(B, {
      ## simulate max stable process with estimated parameters, margins are unit frechet
      x_star <- SpatialExtremes::rmaxstab(nrow(data), coord = locations, cov.mod = Kovmod,  cov11 = cov11_hat,
                                          cov12 = cov12_hat, cov22 = cov22_hat)
      x_star}, simplify = FALSE)
  }
  else {
    if(Kovmod == "brown") {
      nugg_hat <- 0
      range_hat <- parhat[1]
      smooth_hat <- parhat[2]
      }
    else {
      nugg_hat <- parhat[1]
      range_hat <- parhat[2]
      smooth_hat <- parhat[3]
    }

    X_star <- replicate(B, {
      ## simulate max stable process with estimated parameters, margins are unit frechet
      x_star <- SpatialExtremes::rmaxstab(nrow(data), coord = locations, cov.mod = Kovmod,
                                          nugget = nugg_hat,
                                          range = range_hat, smooth = smooth_hat)
      x_star }, simplify = FALSE)

  }

  X_star

}

bootstrap_scalegev  <- function(data, temp.cov, locations,  B = 300, H0 = "ED",
                                ms_models = c("powexp", "gauss", "brown"),
                                sel_crit = "TIC"){

  if(!(H0 %in% c("ED", "LS", "both"))) {
    stop("H0 must be one of 'ED', 'LS' or 'both'.")
  }

  d <- ncol(data)
  nobs <- nrow(data)


  # generate bootstrap samples with unit frechet margins and dependence structure estimated from data

  X_star <- generate_bootsamp_unitfrech(data = data, temp.cov = temp.cov, locations = locations,
                                        ms_models = ms_models, B = B, sel_crit = sel_crit)

  ## estimate parameters of GEV margins under homogeneity assumption
  pars_h0_hom <-  fit_spat_scalegev(data = data, temp.cov = temp.cov, hom = TRUE, returnRatios = FALSE)$mle
  pars_h0_hom[2, ] <- pars_h0_hom[1, ]/pars_h0_hom[2, ]   ## to obtain values for sigma
  pars_h0_hom[4, ] <- pars_h0_hom[4, ]*pars_h0_hom[1, ]   ## to obtain values for alpha
  rownames(pars_h0_hom) <- c("mu0", "sigma0", "gamma0", "alpha0")

  ## estimate GEV parameters on pooled sample
  pars_h0_ed <- fit_gevexpsc(data = as.vector(data),
                             temp.cov = rep(temp.cov, d))$mle
  pars_h0_ed <- matrix(rep(pars_h0_ed, d), ncol = d)
  rownames(pars_h0_ed) <- c("mu0", "sigma0", "gamma0", "alpha0")



  X_star_ed <- purrr::map(X_star, ~ {

    ## transform margins to GEV-distributions with estimated parameters satisfying H0
    for (i in 1:d){
      loc <- pars_h0_ed[1,i]*exp(pars_h0_ed[4,i]/pars_h0_ed[1,i]*temp.cov)
      scale <- pars_h0_ed[2,i]*exp(pars_h0_ed[4,i]/pars_h0_ed[1,i]*temp.cov)

      .x[,i] <-
        SpatialExtremes::frech2gev(.x[,i], loc = loc,
                                   scale = scale,
                                   shape = pars_h0_ed[3,i])
    }

    #  x_star[NAMatrix] <- NA
    .x
  })

  X_star_if <- purrr::map(X_star, ~ {

    ## transform margins to GEV-distributions with estimated parameters satisfying H0
    for (i in 1:d){
      loc <- pars_h0_hom[1,i]*exp(pars_h0_hom[4,i]/pars_h0_hom[1,i]*temp.cov)
      scale <- pars_h0_hom[2,i]*exp(pars_h0_hom[4,i]/pars_h0_hom[1,i]*temp.cov)

      .x[,i] <-
        SpatialExtremes::frech2gev(.x[,i], loc = loc,
                                   scale = scale,
                                   shape = pars_h0_hom[3,i])
    }

    #  x_star[NAMatrix] <- NA
    .x
  })


  # Calculate test statistics for generated bootstrap data
  results_sim_ed <- purrr::map_dfr(X_star_ed,  ~ compute_teststat(.x, temp.cov = temp.cov,
                                                                  equal_distr = TRUE, EDhyp = EDhyp))
  results_sim_if <- purrr::map_dfr(X_star_if,  ~ compute_teststat(.x, temp.cov = temp.cov,
                                                                  equal_distr = FALSE))
  results_ed <- compute_teststat(data = data, temp.cov = temp.cov, equal_distr = TRUE,
                                 EDhyp = EDhyp)
  results_if <- compute_teststat(data = data, temp.cov = temp.cov, equal_distr = FALSE)

  tmp_ed <- rank(c(results_ed$p, results_sim_ed$p), na.last = NA)
  p_corr_ed <- as.vector(tmp_ed[1]/length(tmp_ed))

  tmp_if <- rank(c(results_if$p, results_sim_if$p), na.last = NA)
  p_corr_if <- as.vector(tmp_if[1]/length(tmp_if))

  as_tibble(results_ed)  %>%
    mutate(p_corr = p_corr_ed, bootstrap.B = B, Model_Sel = Kovmod,
           ms_pars = list(parhat), boot_teststat = list(results_sim_ed),
           par_h0 = list(pars_h0_ed), H0 = "Equal Distr.") %>%
    bind_rows(as_tibble(results_if )  %>%
                mutate(p_corr = p_corr_if, bootstrap.B = B, Model_Sel = Kovmod,
                       ms_pars = list(parhat),
                       boot_teststat = list(results_sim_if),
                       par_h0 = list(pars_h0_hom), H0 = "IF"))

}

