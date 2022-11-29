
#' Functionals to express Homogeneity/Equal Distribution hypotheses
#'
#' Functionals \eqn{h} for which, under the null-hypothesis, we have \eqn{h(\theta)  = 0}.
#'
#' @param par.mat  A \eqn{4\times d} matrix containing the estimated parameter vectors for the \eqn{d} stations.
#' @param equal_distr Logical; if `TRUE` (default), the hypothesis is the one of equal distributions across the different
#' stations. If `FALSE`, the hypothesis is that the local scaling model holds across the region.
#'
#' @return A vector of length \eqn{4(d-1)} (equal distribution hypothesis) or \eqn{3(d-1)} (local scaling hypothesis).
#' @export
#'
#' @examples
#' A <- matrix(1:8, ncol = 2)
#' htheta(A)
htheta <- function(par.mat, equal_distr  = TRUE){

  d <- ncol(par.mat)
  if(is.null(d)) {stop("Dimension of parameter matrix must be at least 2.")}

  if(!equal_distr) {
    par.mat <- matrix(c(par.mat[1, ]/par.mat[2,], par.mat[3,], par.mat[4,]/par.mat[1,]),
                      ncol = d,
                      byrow = TRUE)
  }

  as.vector(t(par.mat[ , 1:(d-1)])) -  as.vector(t(par.mat[ , 2:d]))

}

grad_h_delta <- function(par.vec) {

  c(1/par.vec[2], - par.vec[1]/par.vec[2]^2, 0, 0)

}
grad_h_eta <- function(par.vec) {

  c(-par.vec[4]/par.vec[1]^2, 0, 0,  1/par.vec[1])


}

grad_h <- function(par.mat, equal_distr = TRUE) {
  # par.mat must contain parameters mu0, sigma0, gamma0, alpha0 columnwise for the stations
  d <- ncol(par.mat)
  if(is.null(d)) {stop("Dimension of parameter matrix must be at least 2.")}

  if(!equal_distr) {
    abl.delta <- apply(par.mat, 2, grad_h_delta)

    abl.disp <- purrr::map(1:(d-1), ~ c(rep(0, 4*(.x-1)),
                                        abl.delta[, .x], - abl.delta[, .x+1],
                                        rep(0, (4*(d - .x +1) - 8 )) ))


    abl.gam <- purrr::map(1:(d-1), ~ c(rep(0, 4*(.x-1)),
                                       c(0, 0, 1, 0, 0, 0, -1, 0),
                                       rep(0, (4*(d - .x +1) - 8 )) ))

    abl.eta <- apply(par.mat, 2, grad_h_eta)

    abl.eta <- purrr::map(1:(d-1), ~ c(rep(0, 4*(.x-1)),
                                       abl.eta[, .x], - abl.eta[, .x+1],
                                       rep(0, (4*(d - .x +1) - 8 )) ))

    matrix(c(unlist(abl.disp), unlist(abl.gam), unlist(abl.eta) ),
           ncol = 4*d, byrow = TRUE)
  }
  else {

    abl.mu <- purrr::map(1:(d-1), ~ c(rep(0, 4*(.x-1)),
                                      c(1, 0, 0, 0, -1, 0, 0, 0),
                                      rep(0, (4*(d - .x +1) - 8 )) ))

    abl.sigma <- purrr::map(1:(d-1), ~ c(rep(0, 4*(.x-1)),
                                         c(0, 1, 0, 0, 0, -1, 0, 0),
                                         rep(0, (4*(d - .x +1) - 8 )) ))


    abl.gam <- purrr::map(1:(d-1), ~ c(rep(0, 4*(.x-1)),
                                       c(0, 0, 1, 0, 0, 0, -1, 0),
                                       rep(0, (4*(d - .x +1) - 8 )) ))

    abl.alpha <- purrr::map(1:(d-1), ~ c(rep(0, 4*(.x-1)),
                                         c(0, 0, 0, 1, 0, 0, 0, -1),
                                         rep(0, (4*(d - .x +1) - 8 )) ))
    matrix(c(unlist(abl.mu), unlist(abl.sigma), unlist(abl.gam), unlist(abl.alpha)),
           ncol = 4*d, byrow = TRUE)


  }
}


#' Test statistic for testing the hypothesis of Equal Distribution or localing scaling model
#'
#' @param theta Matrix containing estimated parameter values. Each column corresponds to one region.
#' @param covmat Estimated variance-covariance matrix of the estimated parameters.
#' @param n Record length, only needs to be provided if covariance matrix is not scaled by n already.
#' @param H0 Either 'ED' for 'Equal distribution, i.e. the hypothesis is the one of equal distributions across
#' the different stations, or 'LS' for 'Local Scaling', i.e. the hypothesis is that the local scaling model holds across the region.
#' @param covmat_scaled Logical; whether covariance matrix is scaled by record length.
#'
#' @return The numeric value of the test statistic.
#' @export
#'
#' @examples
#' # generate 3 dimensional data:
#' xx <- matrix(rep(exp((1:100)/100), 3)*evd::rgev(300), ncol = 3)
#'
#' # computer scale-GEV parameters for each station and corresponding covariance matrix
#' mlest <- fit_spat_scalegev(data = xx, temp.cov = (1:100)/100)
#'
#' # compute test statistic
#' teststat(mlest$mle, covmat = mlest$cov.mat, n = 100, H0 = "ED")

teststat <- function(theta, covmat, H0 = "ED", covmat_scaled = TRUE, n = NULL) {
  if(anyNA(theta)) {return(NA)}

  if(!covmat_scaled & is.null(n)) { stop("Please provide the record length.")}

  equal_distr <- ifelse(H0 == "ED", TRUE, FALSE)

  hth <- htheta(theta, equal_distr = equal_distr)
  grh <- grad_h(theta,  equal_distr = equal_distr)

  scale.factor <- ifelse(covmat_scaled, 1, n)
  tn <- tryCatch(scale.factor*as.numeric(hth %*%  solve(grh %*% covmat %*% t(grh)) %*% hth),
           error = function(egal) {NA})
  if(is.na(tn)) { warning(" Covariance matrix of h(theta) is singular.")}
  tn

}

# compute tibble containing the value of the test statistic and p-value from the data
#' Compute a data frame containing the value of the test statistic and p-value from the data
#' @param data Matrix or data frame of observations. Each column corresponds to one station/region.
#' @param temp.cov Temporal covariate with same record length as data.
#' @param H0 Either 'ED' for 'Equal distribution, i.e. the hypothesis is the one of equal distributions across
#' the different stations, or 'LS' for 'Local Scaling', i.e. the hypothesis is that the local scaling model holds across the region.
#' @param varmeth Method for estimation of variance-covariance matrix. Can be either `chain` (the default) for an estimator based
#' on the multivariate chain rule, or `basic` for a very simplistic but faster method. Passed to
#' \link[findpoolreg]{fit_spat_scalegev}.
#' @param start_vals Optional: matrix containing start values for the ML optimisation.
#'
#' @return A data frame with columns
#' * teststat: observed value of the test statistic
#' * p: p-value obtained from asymptotic \eqn{\chi^2} distribution (not reliable)
#' @export
#'
#' @examples
#'  # generate 3 dimensional data:
#' xx <- matrix(rep(exp((1:100)/100), 3)*evd::rgev(300), ncol = 3)
#' compute_teststat(data = xx, temp.cov = (1:100)/100, H0 = "ED")
compute_teststat <- function(data, temp.cov, H0 = "ED", varmeth = "chain", start_vals = NULL) {

  d <- ncol(data)
  equal_distr <- ifelse(H0 == "ED", TRUE, FALSE)

  if(is.null(d)) {stop("Data must have dimension at least 2.")}
  mlest <- tryCatch(fit_spat_scalegev(data, temp.cov, varmeth = varmeth, start_vals = start_vals),
                    error = function(egal) list(mle = NA, cov.mat = NA))

  val <- teststat(theta = mlest$mle, covmat = mlest$cov.mat, n = nrow(data),
                  H0 = H0, covmat_scaled = TRUE)

  dftest <- ifelse(equal_distr, 4*(d-1), 3*(d-1))
  data.frame(teststat = val, p = stats::pchisq(val, df = dftest,
                                                lower.tail = FALSE))
}
