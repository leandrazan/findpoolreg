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
#' @return Returns a list of  3.
#' * `X_star` contains the `B` bootstrap samples, each being a matrix of same dimension as input data.
#' * `covmod` is the covariance model of the chosen max-stable process
#' * `parhat` the parameters of the chosen max-stable process
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
#' plot.ts(bootsamps$X_star[[1]])
#' plot.ts(bootsamps$X_star[[2]])
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
    fit.tmp <-  tryCatch(
      SpatialExtremes::fitmaxstab(unitfrech, coord = locations, ms.mod, warn = warn_msfit),
        error = function(a) NA)
    fitted_ms <- if(!is.na(fit.tmp)) {
      append(fitted_ms, list(fit.tmp))
      }
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

  list(X_star = X_star, covmod = Kovmod, parhat = list(parhat))

}


#' Bootstrap p-values for test-statistic
#'
#' @param data A matrix of observations. Each column corresponds to one station.
#' @param temp.cov A numeric vector containing values of the temporal covariate. Must be of length `nrow(data)`.
#' @param locations A matrix containing the coordinates (longitude and latitude) of the stations.
#' The i-th row corresponds to the data contained in the i-th column of the data matrix.
#' @param B The number of bootstrap replicates.
#' @param H0 Character code for the null-hypothesis: `ED` for the one of equal distribution;
#' `LS` for the one of a local scaling model.
#' @param ms_models A vector containing the names of the max-stable models to fit. Must be a subset of
#' `gauss, brown, powexp, cauchy, whitmat, bessel`.
#' @param sel_crit The selection criterion based on which the max-stable model is selected. Must be one of
#' `TIC, AIC`.
#' @param varmeth Method for estimating the variance-covariance matrix of the stationwise ML estimators, passed to
#' \code{\link[findpoolreg]{fit_spat_scalegev}}.
#'  Can be either `chain` (the default) for an estimator based
#' on the multivariate chain rule, or `basic` for a very simplistic but faster method.
#'
#' @return A tibble with the following columns:
#' * `teststat` The value of the test statistic on the observed data
#' * `p` The asymptotic p value of the test statistic
#' * `p_boot` The bootstrapped p-value.
#' * `bootstrap.B` The number of bootstrap repetitions that were performed.
#' * `Model_Sel` The max-stable process model that was selected.
#' * `ms_pars` Parameters of the selected max-stable process model.
#' * `boot_teststat` Data frame containing the test statistic and corresponding p-value for each bootstrap replicate.
#' * `par_h0` The parameter estimates under H0.
#' * `H0` Character code of the null hypothesis.
#' @export
#'
#' @examples
#' \dontrun{
#' # generate some data
#' cvrt <- (1:100)/100
#' coords <- matrix(20*abs(stats::rnorm(4*2)), ncol = 2)
#' x <- generateData(seed = 2, n = 100, temp.cov = cvrt, d= 4, locations = coords)
#'
#' bootres <- bootstrap_scalegev(data = x, temp.cov = cvrt, locations = coords, varmeth = "chain", B = 200)
#' bootres
#' }
bootstrap_scalegev  <- function(data, temp.cov, locations,  B = 300, H0 = "ED",
                                ms_models = c("powexp", "gauss", "brown"),
                                sel_crit = "TIC", varmeth = "chain"){

  if(!(H0 %in% c("ED", "LS", "both"))) {
    stop("H0 must be one of 'ED', 'LS' or 'both'.")
  }

  d <- ncol(data)
  nobs <- nrow(data)


  # generate bootstrap samples with unit frechet margins and dependence structure estimated from data

  X_star <- generate_bootsamp_unitfrech(data = data, temp.cov = temp.cov, locations = locations,
                                        ms_models = ms_models, B = B, sel_crit = sel_crit)
  covmod <- X_star$covmod
  parhat_ms <- X_star$parhat
  X_star <- X_star$X_star

  if(H0 %in% c("LS", "both")) {

    ## estimate parameters of GEV margins under homogeneity assumption
    pars_h0_hom <-  fit_spat_scalegev(data = data, temp.cov = temp.cov, hom = TRUE, returnRatios = FALSE)$mle
    pars_h0_hom[2, ] <- pars_h0_hom[1, ]/pars_h0_hom[2, ]   ## to obtain values for sigma
    pars_h0_hom[4, ] <- pars_h0_hom[4, ]*pars_h0_hom[1, ]   ## to obtain values for alpha
    rownames(pars_h0_hom) <- c("mu0", "sigma0", "gamma0", "alpha0")

    X_star_if <- purrr::map(X_star, ~ {

      ## transform margins to GEV-distributions with estimated parameters satisfying H0
      for (i in 1:d){
        expo <- exp(pars_h0_hom[4,i]/pars_h0_hom[1,i]*temp.cov)
        loc <- pars_h0_hom[1,i]*expo
        scale <- pars_h0_hom[2,i]*expo

        .x[,i] <-
          SpatialExtremes::frech2gev(.x[,i], loc = loc,
                                     scale = scale,
                                     shape = pars_h0_hom[3,i])
      }

      #  x_star[NAMatrix] <- NA
      .x
    })

    results_sim_if <- purrr::map_dfr(X_star_if,  ~ compute_teststat(.x, temp.cov = temp.cov,
                                                                    H0 = "LS", varmeth = varmeth))
    results_if <- compute_teststat(data = data, temp.cov = temp.cov, H0 = "LS", varmeth = varmeth)

    tmp_if <- rank(c(results_if$p, results_sim_if$p), na.last = NA)
    p_corr_if <- as.vector(tmp_if[1]/length(tmp_if))

    tibif <- dplyr::as_tibble(results_if )  %>%
                dplyr::mutate(p_boot = p_corr_if, bootstrap.B = B, Model_Sel = covmod,
                       ms_pars = parhat_ms,
                       boot_teststat = list(results_sim_if),
                       par_h0 = list(pars_h0_hom), H0 = "LS")
  }

  if(H0 %in% c("ED", "both")) {
    ## estimate GEV parameters on pooled sample
    pars_h0_ed <- fit_scalegev(data = data,
                               temp.cov = temp.cov, hessian = FALSE)$mle
    pars_h0_ed <- matrix(rep(pars_h0_ed, d), ncol = d)
    rownames(pars_h0_ed) <- c("mu0", "sigma0", "gamma0", "alpha0")

    X_star_ed <- purrr::map(X_star, ~ {

      ## transform margins to GEV-distributions with estimated parameters satisfying H0
      for (i in 1:d){
        expo <- exp(pars_h0_ed[4,i]/pars_h0_ed[1,i]*temp.cov)
        loc <- pars_h0_ed[1,i]*expo
        scale <- pars_h0_ed[2,i]*expo

        .x[,i] <-
          SpatialExtremes::frech2gev(.x[,i], loc = loc,
                                     scale = scale,
                                     shape = pars_h0_ed[3,i])
      }

      #  x_star[NAMatrix] <- NA
      .x
    })

    # Calculate test statistics for generated bootstrap data
    results_sim_ed <- purrr::map_dfr(X_star_ed,  ~ compute_teststat(.x, temp.cov = temp.cov,
                                                                    H0 = "ED", varmeth = varmeth))

    results_ed <- compute_teststat(data = data, temp.cov = temp.cov, H0 = "ED", varmeth = varmeth)

    tmp_ed <- rank(c(results_ed$p, results_sim_ed$p), na.last = NA)
    p_corr_ed <- as.vector(tmp_ed[1]/length(tmp_ed))

   tibed <-  dplyr::as_tibble(results_ed)  %>%
     dplyr::mutate(p_boot = p_corr_ed, bootstrap.B = B, Model_Sel = covmod,
                   ms_pars = parhat_ms,
                   boot_teststat = list(results_sim_ed),
                   par_h0 = list(pars_h0_ed), H0 = "ED")

  }


  if(H0 == "both") {
    return( tibed %>% dplyr::bind_rows(tibif))
  } else if ( H0 == "ED") {
    return(tibed)
  } else {
    return(tibif)
  }

}




#' Bootstrap p-values of test-statistic on several subsets
#'
#' @param data A matrix of observations. Each column corresponds to one station.
#' Columns should be named or numbered to identify subsets.
#' @param subsets A list containing subsets of columns on which to perform the bootstrap test.
#' @param return_boots Logical; whether to return bootstrapped values of test statistic and p-value.
#' @param adj_pvals Logical; whether to return column containing Holmes-adjusted p-values.
#' @inheritParams bootstrap_scalegev
#'
#' @return A tibble with the following columns:
#' * `teststat` The value of the test statistic on the observed data
#' * `p` The asymptotic p value of the test statistic
#' * `p_boot` The bootstrapped p-value.
#' * `bootstrap.B` The number of bootstrap repetitions that were performed.
#' * `Model_Sel` The max-stable process model that was selected.
#' * `ms_pars` Parameters of the selected max-stable process model.
#' * `boot_teststat` (if `return_boots = TRUE`) Data frame containing the test statistic and corresponding p-value for each bootstrap replicate.
#' * `par_h0` The parameter estimates under H0.
#' * `H0` Character code of the null hypothesis.
#' * `sbst` The subset on which the hypothesis was tested.
#' * `adj_p` (if `adj_pvals = TRUE`) The Holmes-adjuted p-values.
#' @export
#'
#' @examples
#'
#' \dontrun{
#' # generate 6-dimensional data where two stations deviate from the rest
#' cvrt <- (1:100)/100
#' coords <- matrix(20*abs(stats::rnorm(6*2)), ncol = 2)
#' x <- generateData(seed = 2, n = 100, temp.cov = cvrt, d= 6, locations = coords,
#' loc = c(rep(10,4), c(8, 8)), scale = c(rep(2,4), c(1.5, 1.5)), shape = rep(0.1, 6), alpha = rep(2,6))
#'
#' colnames(x) <- 1:6
#' sbsts <- list( c(1,2), c(1,3), c(1,4), c(1,5), c(1,6))
#' bootres <- bootstrap_scalegev_subsets(data = x, temp.cov = cvrt, locations = coords, varmeth = "chain", B = 200, subsets = sbsts,
#' adj_pvals = TRUE)
#' bootres %>% tidyr::unnest(cols = sbst) # H0 ist rejected at 5%-level for the deviating stations 5 and 6
#' }
bootstrap_scalegev_subsets  <- function(data, temp.cov, locations,  B = 300, H0 = "ED",
                                ms_models = c("powexp", "gauss", "brown"),
                                sel_crit = "TIC", varmeth = "chain", subsets, return_boots = FALSE, adj_pvals = FALSE){

  if(!(H0 %in% c("ED", "LS", "both"))) {
    stop("H0 must be one of 'ED', 'LS' or 'both'.")
  }
  if(missing(subsets)) { stop("Please provide the subsets on which to perform the test.")}

  d <- ncol(data)
  nobs <- nrow(data)

  if(is.null(colnames(data))) { colnames(data) <- 1:d}

  # generate bootstrap samples with unit frechet margins and dependence structure estimated from data

  X_star <- generate_bootsamp_unitfrech(data = data, temp.cov = temp.cov, locations = locations,
                                        ms_models = ms_models, B = B, sel_crit = sel_crit)

  covmod <- X_star$covmod
  parhat_ms <- X_star$parhat
  X_star <- X_star$X_star

  res <- dplyr::tibble()
  # compute bootstrap test statistics  for each subset

  for( j in 1:length(subsets)) {


    j.cols <-  subsets[[j]]

    # observations belonging to current subset
    data.j <- data[, colnames(data) %in% j.cols]

    ## estimate GEV parameters on pooled sample consisting of current data
    pars_h0_ed <- fit_scalegev(data = data.j,
                               temp.cov = temp.cov)$mle

    # get columns of bootstrap samples corresponding to current subset
    X_star_ed <- purrr::map(X_star, ~ { .x[ , (colnames(data) %in% j.cols)]})

    # H0 parameter estimates of current subset
    expo <- exp(pars_h0_ed[4]/pars_h0_ed[1]*temp.cov)
    loc <- pars_h0_ed[1]*expo
    scale <- pars_h0_ed[2]*expo

    # transform margins of bootstrap samples to GEV margins satisfying H0
    X_star_ed <- purrr::map(X_star_ed, ~ {

      ## transform margins to GEV-distributions with estimated parameters satisfying H0
      for (i in  1:length(j.cols)){

        .x[, i] <-
          SpatialExtremes::frech2gev(.x[, i], loc = loc,
                                     scale = scale,
                                     shape = pars_h0_ed[3])
      }

      #  x_star[NAMatrix] <- NA
      .x
    })

    # ## estimate parameters of GEV margins under homogeneity assumption
    # pars_h0_hom <-  fit_spat_expsc(data = data.j, temp.cov = temp.cov, hom = TRUE)$mle
    # pars_h0_hom[2, ] <- pars_h0_hom[1, ]/pars_h0_hom[2, ]   ## to obtain values for sigma
    # pars_h0_hom[4, ] <- pars_h0_hom[4, ]*pars_h0_hom[1, ]   ## to obtain values for alpha
    # rownames(pars_h0_hom) <- c("mu0", "sigma0", "gamma0", "alpha0")


    # X_star_if <- purrr::map(X_star, ~ {
    #
    #   ## transform margins to GEV-distributions with estimated parameters satisfying H0
    #   for (i in 1:length(j.cols)){
    #     loc <- pars_h0_hom[1,i]*exp(pars_h0_hom[4,i]/pars_h0_hom[1,i]*temp.cov)
    #     scale <- pars_h0_hom[2,i]*exp(pars_h0_hom[4,i]/pars_h0_hom[1,i]*temp.cov)
    #
    #     .x[, (colnames(data) == j.cols[i])] <-
    #       SpatialExtremes::frech2gev(.x[,(colnames(data) == j.cols[i])], loc = loc,
    #                                  scale = scale,
    #                                  shape = pars_h0_hom[3,i])
    #   }
    #
    #   #  x_star[NAMatrix] <- NA
    #   .x[, (colnames(data) %in% j.cols)]
    # })


    # Calculate test statistics for generated bootstrap data
    results_sim_ed <- purrr::map_dfr(X_star_ed,  ~ compute_teststat(.x, temp.cov = temp.cov,
                                                                   H0 = "ED", varmeth = varmeth))

    # results_sim_if <- purrr::map_dfr(X_star_if,  ~ compute_teststat(.x, temp.cov = temp.cov,
    #                                                                 equal_distr = FALSE))
    results_ed <- compute_teststat(data = data.j, temp.cov = temp.cov,
                                   H0 = "ED", varmeth = varmeth)
    #results_if <- compute_teststat(data = data.j, temp.cov = temp.cov, equal_distr = FALSE)

    tmp_ed <- rank(c(results_ed$p, results_sim_ed$p), na.last = NA)
    p_corr_ed <- as.vector((tmp_ed[1]-1)/length(tmp_ed))

    # tmp_if <- rank(c(results_if$p, results_sim_if$p), na.last = NA)
    # p_corr_if <- as.vector((tmp_if[1]-1)/length(tmp_if))

    if(return_boots) {
      res <- res %>% dplyr::bind_rows(
        dplyr::as_tibble(results_ed)  %>%
          dplyr::mutate(p_boot = p_corr_ed, bootstrap.B = B, Model_Sel = covmod,
                 ms_pars = parhat_ms, boot_teststat = list(results_sim_ed),
                 par_h0 = list(pars_h0_ed), H0 = "ED", sbst  = list(data.frame(t(subsets[[j]])))))
        # %>%
        #   bind_rows(as_tibble(results_if )  %>%
        #               mutate(p_corr = p_corr_if, bootstrap.B = B, Model_Sel = Kovmod,
        #                      ms_pars = list(parhat),
        #                      boot_teststat = list(results_sim_if),
        #                      par_h0 = list(pars_h0_hom), H0 = "IF"))  %>% mutate( sbst  = list(data.frame(t(subsets[[j]])))))
    }
    else {
      res <- res %>% dplyr::bind_rows(
        dplyr::as_tibble(results_ed)  %>%
          dplyr::mutate(p_boot = p_corr_ed, bootstrap.B = B, Model_Sel = covmod,
                 ms_pars = parhat_ms,
                 par_h0 = list(pars_h0_ed), H0 ="ED", sbst  = list(data.frame(t(subsets[[j]])))))
      # %>%
      #     bind_rows(as_tibble(results_if )  %>%
      #                 mutate(p_corr = p_corr_if, bootstrap.B = B, Model_Sel = Kovmod,
      #                        ms_pars = list(parhat),
      #                        par_h0 = list(pars_h0_hom), H0 = "IF"))  %>% mutate(sbst  = list(data.frame(t(subsets[[j]])))))
    }
  }

  if(adj_pvals) {
    adjp <- adjust_pvals(res$p_boot)
    res <- res %>% dplyr::arrange(p_boot) %>% dplyr::mutate(adj_p = adjp)
  }

  res

}




#' Holms-correction of p-values for multiple testing
#'
#' @param pvals  Vector of p-values
#'
#' @return Vector of adjusted p-values
#' @export
#'
#' @examples
#' adjust_pvals( c(0.001, 0.0002, 0.3, 0.55) )
adjust_pvals <- function(pvals) {

  n.na <- sum(is.na(pvals))
  pvals <- pvals[!is.na(pvals)]
  nhyp <- length(pvals)

  pvals_sort <- sort(pvals, na.last = TRUE)

  adjp <- (nhyp:1)*pvals_sort

  for( j in 2:nhyp) {
    adjp[j] <- max(adjp[j-1], adjp[j])
  }

  adjp[adjp > 1] <-  1

  c(adjp, rep(NA, n.na))

}


