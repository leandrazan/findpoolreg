
# computes \dot \ell_{(0, 1, \gamma)}(( x- mut)/sigmat)
#' Computes three dimensional score function of data with margins transformed to standard GEV
#'
#' @param x Numeric vector of observations from the scale-GEV-model
#' @param theta Parameter vector of the scale-GEV-model
#' @param temp.cov Values of temporal covariate; numeric vector of same length as x
#' @param rel_par Logical; whether `theta` contains the plain parameters \eqn{\mu, \sigma, \gamma, \alpha},
#' or relative values \eqn{\mu, \mu/\sigma, \gamma, \alpha/\mu}.
#'
#' @return A matrix of dimension \eqn{3\times} `length(x)` containing the values of
#' \deqn{ \dot\ell_{(0,1, \gamma)} ( \frac{ x - \mu(gmst(t))}{\sigma(gmst(t))} ).}
#' @export
#'
#' @examples
#'  xx <- exp((1:100)/100)*evd::rgev(100)
#'  score_standard_univ(x = xx, theta = c(10, 2, 0, 1), temp.cov = (1:100)/100, rel_par = FALSE)
score_standard_univ <- function(x, theta, temp.cov = NULL, rel_par = FALSE) {

  if(is.null(temp.cov)) {stop("Temporal covariate must be provided.")}

  # change to mu sigma, gamma, alpha notation
  if(rel_par) {
    mu0 <- theta[1]
    sigma0 <- theta[1]/theta[2]
    xi <- theta[3]
    alpha0 <- theta[4]*theta[1]
  }

  else {
    mu0 <- theta[1]
    sigma0 <- theta[2]
    xi <- theta[3]
    alpha0 <- theta[4]
  }


  mut <- mu0*exp(alpha0*temp.cov/mu0)
  sigmat <-  sigma0*exp(alpha0*temp.cov/mu0)

  zt <- (x - mut)/sigmat
  in_supp <- which(1+ xi*zt > 0)
  zt[!in_supp] <- NA

  if(abs(xi) < 1e-08) {
    ut <- exp(-zt)
  }
  else {
    ut <- (1+xi *zt)^(-1/xi)
  }

  scoreloc <- (xi + 1 - ut)/((1+xi*zt))
  scorescale <- ((1-ut)*zt -1)/((1+xi*zt))

  if(abs(xi) < 1e-08){
    scoreshape <- (1-ut)*zt^2/2- zt
  }
  else{
    scoreshape <- (1-ut)*1/xi*(1/xi*log(1+xi*zt)-zt/(1+xi*zt)) - zt/(1+xi*zt)
  }

  scoremat <- matrix(c(scoreloc, scorescale, scoreshape),
                     nrow = 3, byrow = TRUE)

  return(scoremat)

}


#' Computes cross-covariance matrix of score function applied to columns of
#' data transformed to standard GEV distribution
#'
#' The empirical covariance matrices of
#' \deqn{ (\dot\ell_{(0, 1, \hat\gamma_j)}( (x_t - \mu_j(gmst(t)) )/\sigma_j(gmst(t))) )_t \text{ and }
#' (\dot\ell_{(0, 1, \hat\gamma_k)}( (x_t - \mu_k(gmst(t)) )/\sigma_k(gmst(t))) )_t}
#'
#'  are needed for each combination of \eqn{j, k} in the estimation of the ML estimator's covariance matrix.
#'
#'
#' @param data Matrix of observations.
#' @param parmat Matrix of parameters of the scale-GEV-model. Each column corresponds to one station.
#' @param temp.cov  Values of temporal covariate; numeric vector of length as `nrow(data)`.
#' @param rel_par Logical; whether `theta` contains the plain parameters \eqn{\mu, \sigma, \gamma, \alpha},
#' or relative values \eqn{\mu, \mu/\sigma, \gamma, \alpha/\mu}.
#'
#' @return A \eqn{d\times d} matrix.
#' @export
#'
#' @examples
#' # generate data
#' cvrt <- (1:50)/50
#' d <- 4
#' coords <- matrix(20*abs(stats::rnorm(d*2)), ncol = 2)
#' x <- generateData(seed = 1, n = 50, temp.cov = cvrt, d= d, locations = coords)
#' mle <- fit_spat_scalegev(x, temp.cov =cvrt)
#' estimate_gammas(x, parmat = mle$mle, temp.cov =cvrt, rel_par = FALSE)
estimate_gammas <- function(data, parmat, temp.cov = NULL, rel_par = FALSE) {
  n <- nrow(data)
  d <- ncol(data)
  scores <- array(dim = c(3, n, d))
  Gammas <- array(dim = c(3*d,3*d))
  for( j in 1:d ) {
    scores[ , , j] <- score_standard_univ(data[, j], theta = parmat[, j], temp.cov = temp.cov, rel_par = rel_par)
  }
  for (j in 1:d) {
    for (k in j:d) {
      Gammas[ ((j -1)*3 +1) : (j*3), ((k -1)*3 +1) : (k*3) ] <- cov(t(scores[, , j]), t(scores[ , , k]))
      Gammas[ ((k -1)*3 +1) : (k*3), ((j -1)*3 +1) : (j*3) ] <- t(Gammas[ ((j -1)*3 +1) : (j*3), ((k -1)*3 +1) : (k*3) ] )
    }
  }
  Gammas
}


# compute entries of matrices B_c(mu, sigma, alpha)
Bmat <- function(par, temp.cov) {

  n <- length(temp.cov)

  mu <- par[1]
  sigma <- par[2]
  alpha <- par[4]

  expo <- exp(alpha/mu * temp.cov)

  b11 <- ( 1- alpha*temp.cov/mu)*expo
  b12 <- - sigma*alpha*temp.cov/mu^2*expo
  b13 <- rep(0, n)

  b21 <- rep(0, n)
  b22 <- expo
  b23 <- rep(0, n)

  b31 <- b32 <- rep(0, n)
  b33 <- rep(1, n)

  b41 <- temp.cov* expo
  b42 <- sigma*temp.cov/mu * expo
  b43 <- rep(0, n)

 assemble_mat(b11, b12, b13, b21, b22, b23, b31, b32, b33, b41, b42, b43)

}

# assemble the matrices B_c(mu, sigma, alpha) for each value of c
assemble_mat <- function(b11, b12, b13, b21, b22, b23, b31, b32, b33, b41, b42, b43) {
  Bs <- list()

  for( j in 1:length(b11)) {

   Bs[[j]] <- matrix( c(b11[j], b12[j], b13[j],
                         b21[j], b22[j], b23[j],
                         b31[j], b32[j], b33[j],
                        b41[j], b42[j], b43[j]),
                      ncol = 3, byrow = TRUE)
  }
  Bs

}

# compute the matrices T_sigma(c)^(-1) for each value of c
Tsigma_inv <- function(par, temp.cov) {

  mu <- par[1]
  sigma <- par[2]
  alpha <- par[4]
  sigmat <- sigma*exp(alpha/mu*temp.cov)

  purrr::map(sigmat, ~ diag(c(1/.x, 1/.x, 1) ))

}

# compute Sigma_{jk}

#' Estimate covariance of ML estimation applied to two stations
#'
#' @param par.j Estimated parameter vector for station j
#' @param par.k Estimated parameter vector for station k
#' @param Gammajk Estimated covariance matrix of score function applied to data standardised to GEV margins,
#' as computed by \code{\link[findpoolreg]{estimate_gammas}}
#' @param temp.cov Values of temporal covariate
#' @param Jinv.j Inverse of hessian returned by calling \code{\link[findpoolreg]{fit_scalegev}} on data from station j
#' @param Jinv.k Inverse of hessian returned by calling \code{\link[findpoolreg]{fit_scalegev}} on data from station k
#'
#' @return A 4x4 matrix containing an estimation of
#' \deqn{ \Cov( \hat \vartheta_j, \hat\vartheta_k).}
#' @export
#'
#' @examples
#'
#' # generate data
#' cvrt <- (1:50)/50
#' d <- 4
#' coords <- matrix(20*abs(stats::rnorm(d*2)), ncol = 2)
#' x <- generateData(seed = 1, n = 50, temp.cov = cvrt, d= d, locations = coords)
#'
#' # compute ML estimate at each station
#' mle1 <- fit_scalegev(x[, 1], temp.cov =cvrt)
#' mle3 <- fit_scalegev(x[, 3], temp.cov =cvrt)
#'
#'
#' # compute covariance between ML estimate of station 1 and ML estimate of station 3
#' Gams <-  estimate_gammas(x[, c(1,3)], parmat = cbind(mle1$mle, mle3$mle), temp.cov =cvrt, rel_par = FALSE)
#' compute_sigmajk(mle1$mle, mle3$mle, Gams[ 1:3, 4:6], cvrt, solve(mle1$hessian), solve(mle3$hessian))

compute_sigmajk <- function(par.j, par.k, Gammajk, temp.cov, Jinv.j , Jinv.k) {

  Bj <- Bmat(par.j, temp.cov = temp.cov)
  Bk <- Bmat(par.k, temp.cov = temp.cov)

  Tj <- Tsigma_inv(par.j, temp.cov = temp.cov)
  Tk <- Tsigma_inv(par.k, temp.cov = temp.cov)

  BmalT.j <- purrr::map2(Bj, Tj, ~ {.x %*% .y})

  BmalT.k <- purrr::map2(Bk, Tk, ~ {.x %*% .y})

  Summands <- purrr::map2(BmalT.j, BmalT.k, ~ {.x %*% Gammajk %*% t(.y)} )

  Summands <- Reduce("+", Summands)/length(temp.cov)

  Jinv.j %*%  Summands %*% Jinv.k

}
