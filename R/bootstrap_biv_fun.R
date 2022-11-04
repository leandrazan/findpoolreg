
generate_bootsamp_unitfrech_bivdens <- function(data, temp.cov, locations,
                                        biv_models = c( "log", "alog", "hr"), B = 500) {

  if(!all(biv_models %in% c( "log", "alog", "hr", "neglog", "aneglog", "bilog",
                            "negbilog", "ct","amix"))) {
    stop("Bivariate extreme value distributions can be \n
   'log', 'alog', 'hr', 'neglog', 'aneglog', 'bilog',
                                                        'negbilog', 'ct','amix').")
  }


  ## transform margins to unit Frechet based on local GEV parameters
  unitfrech <- spat_scalegev2frech(data = data, temp.cov = temp.cov)

  ## Fit max-stable processes using models specified in ms_models (margins are unit Frechet)
  fitted_biv <- list()
  for( biv.mod in biv_models) {
    fit.tmp <-  tryCatch(
     evd::fbvevd(unitfrech, model = biv.mod, loc1 = 1, loc2 = 1, scale1 = 1, scale2 = 1, shape1 = 1, shape2 = 1, std.err = FALSE),
      error = function(a) {list(estimate  =  NA)} )
    fitted_biv <- if(!anyNA(fit.tmp$estimate)) {
      append(fitted_biv, list(fit.tmp))
    }
  }

  ## choose best fit based on selection criterion
  sel.vals <- purrr::map_dbl(fitted_biv, ~ AIC(.x))

  model.sel <- which.min(sel.vals)

  # chosen fitted max-stable model
  fitted <- fitted_biv[[model.sel]]

  # covariance model of chosen bivariate EV model
  biv.mod <- fitted$model

  ## estimated parameters of the bivariate EV model
  parhat <- fitted$estimate

  # generate bootstrap samples
  if(biv.mod %in% c("hr", "log", "neglog")) {

    # sim data
    X_star <- replicate(B, {
      ## simulate max stable process with estimated parameters, margins are unit frechet
      x_star <- evd::rbvevd(nrow(data), dep = parhat, model = biv.mod, mar1 = c(1,1,1), mar2 = c(1,1,1))
      x_star}, simplify = FALSE)

  }
  else if(biv.mod %in% c("alog", "aneglog")) {

    # sim data
    X_star <- replicate(B, {
      ## simulate max stable process with estimated parameters, margins are unit frechet
      x_star <- evd::rbvevd(nrow(data), dep = parhat[3], asy = parhat[1:2], model = biv.mod, mar1 = c(1,1,1), mar2 = c(1,1,1))
      x_star}, simplify = FALSE)
  }
  else  {

    # sim data
    X_star <- replicate(B, {
      ## simulate max stable process with estimated parameters, margins are unit frechet
      x_star <- evd::rbvevd(nrow(data), alpha = parhat[1], beta = parhat[2], model = biv.mod, mar1 = c(1,1,1), mar2 = c(1,1,1))
      x_star}, simplify = FALSE)

  }

  list(X_star = X_star, bivmod = biv.mod, parhat = list(parhat))

}


bootstrap_onepair_bivmod  <- function(data, temp.cov, B = 300, H0 = "ED",
                                        biv_models = c( "log", "alog", "hr"),
                                        varmeth = "chain", return_boots = FALSE){

  if(!(H0 %in% c("ED", "LS", "both"))) {
    stop("H0 must be one of 'ED', 'LS' or 'both'.")
  }

  d <- ncol(data)
  nobs <- nrow(data)


  # generate bootstrap samples with unit frechet margins and dependence structure estimated from data

  X_star <- generate_bootsamp_unitfrech_bivdens(data = data, temp.cov = temp.cov,
                                        biv_models = biv_models, B = B)
  covmod <- X_star$bivmod
  parhat_biv <- X_star$parhat
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

        .x[,i] <- frech2gev(.x[,i], loc = loc,
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
                    biv_pars = parhat_biv,
                    par_h0 = list(pars_h0_hom), H0 = "LS")
    if(return_boots) {tibif <- tibif %>% dplyr::mutate( boot_teststat = list(results_sim_if))}
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

        .x[,i] <- frech2gev(.x[,i], loc = loc,
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
                    biv_pars = parhat_biv,
                    par_h0 = list(pars_h0_ed), H0 = "ED")

    if(return_boots) { tibed <- tibed %>% dplyr::mutate( boot_teststat = list(results_sim_ed))}
  }


  if(H0 == "both") {
    return( tibed %>% dplyr::bind_rows(tibif))
  } else if ( H0 == "ED") {
    return(tibed)
  } else {
    return(tibif)
  }

}




#' Bootstrap p-values of test-statistic on pairs
#'
#' @param data A matrix of observations. Each column corresponds to one station.
#' Columns should be named or numbered to identify subsets.
#' @param temp.cov  A numeric vector containing values of the temporal covariate. Must be of length `nrow(data)`.
#' @param B The number of bootstrap replicates.
#' @param H0 Character code for the null-hypothesis: `ED` for the one of equal distribution;
#' `LS` for the one of a local scaling model.
#' @param biv_models A vector containing the names of the max-stable models to fit. Must be a subset of
#' `gauss, brown, powexp, cauchy, whitmat, bessel`.
#' @param varmeth Method for estimating the variance-covariance matrix of the stationwise ML estimators, passed to
#' \code{\link[findpoolreg]{fit_spat_scalegev}}.
#'  Can be either `chain` (the default) for an estimator based
#' on the multivariate chain rule, or `basic` for a very simplistic but faster method.
#' @param reg_of_int The region of primary interest for the analysis. Can be omitted if `pairs` is provided.
#' @param pairs A list of pairs on which to perform the test.
#' @param return_boots Logical; whether to return bootstrapped test statistics and p values.
#'
#' @return  A tibble with the following columns:
#' * `teststat` The value of the test statistic on the observed data
#' * `p` The asymptotic p value of the test statistic
#' * `p_boot` The bootstrapped p-value.
#' * `bootstrap.B` The number of bootstrap repetitions that were performed.
#' * `Model_Sel` The max-stable process model that was selected.
#' * `biv_pars` Parameters of the selected max-stable process model.
#' * `boot_teststat` (if `return_boots = TRUE`) Data frame containing the test statistic and corresponding p-value for each bootstrap replicate.
#' * `par_h0` The parameter estimates under H0.
#' * `H0` Character code of the null hypothesis.
#' * `X1`, `X2` The labels of pairs on which the hypothesis was tested.
#' * `adj_p` (if `adj_pvals = TRUE`) The Holm-adjusted p-values.
#' @export
#'
#' @examples
#' \dontrun{
#' # generate some data
#' cvrt <- (1:100)/100
#' coords <- matrix(20*abs(stats::rnorm(4*2)), ncol = 2)
#' x <- generateData(seed = 2, n = 100, temp.cov = cvrt, d= 4, locations = coords)
#' bootres <- bootstrap_pairs_bivmod(data = x, temp.cov = cvrt, varmeth = "chain", B = 200, reg_of_int = 3)
#' bootres
#' }
bootstrap_pairs_bivmod <- function(data, temp.cov, B = 300, H0 = "ED",
                                    biv_models = c( "log", "alog", "hr"),
                                    varmeth = "chain", reg_of_int = NULL, pairs = NULL,
                                    return_boots = FALSE) {

  if(is.null(pairs)) {
    if(is.null(reg_of_int)) {
      stop("Please provide either the region of intereset or a list of pairs on which to perform the test.")
    }
    else {
      regs <- colnames(data)
      regs <- regs[!(regs == reg_of_int)]
      pairs <- purrr::map(regs, ~ c(reg_of_int, .x))
    }
  }

  n.pairs <- length(pairs)
  bootres <- tibble::tibble()
  for(j.p in 1:n.pairs) {
    pair.tmp <- pairs[[j.p]]

    data.tmp <- data[, pair.tmp]

    boot.tmp <- bootstrap_onepair_bivmod(data = data.tmp, temp.cov = temp.cov, B = B, H0 = H0,
                                          biv_models = biv_models, return_boots = return_boots)

    bootres <- bootres %>%
      dplyr::bind_rows(boot.tmp %>%
                         dplyr::mutate(X1 = pair.tmp[1], X2 = pair.tmp[2]))

    }

  bootres
}
