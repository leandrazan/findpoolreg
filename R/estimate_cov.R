
# computes \dot \ell_{(0, 1, \gamma)}(( x- mut)/sigmat)

#' Computes three dimensional score function of data with margins transformed to
#' standard GEV
#'
#' @param x Numeric vector of observations from the scale-GEV-model
#' @param theta Parameter vector of the scale-GEV-model
#' @param temp.cov Values of temporal covariate; numeric vector of same length as x
#' @param rel_par Logical; whether `theta` contains the plain parameters \eqn{\mu, \sigma, \gamma, \alpha},
#' or relative values \eqn{\mu, \mu/\sigma, \gamma, \alpha/\mu}.
#'
#' @return A matrix of dimension  `3 x length(x)`.
#' @export
#'
#' @examples
#'  xx <- exp((1:100)/100)*evd::rgev(100)
#'  score_standard_univ(x = xx, theta = c(10, 2, 0, 1),
#'                      temp.cov = (1:100)/100, rel_par = FALSE)
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
#' @param data Matrix of observations.
#' @param parmat Matrix of parameters of the scale-GEV-model. Each column corresponds to one station.
#' @param temp.cov  Values of temporal covariate; numeric vector of length as `nrow(data)`.
#' @param rel_par Logical; whether `theta` contains the plain parameters
#' \eqn{\mu, \sigma, \gamma, \alpha},
#' or relative values \eqn{\mu, \mu/\sigma, \gamma, \alpha/\mu}.
#'
#' @return A dxd matrix.
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
      Gammas[ ((j -1)*3 +1) : (j*3), ((k -1)*3 +1) : (k*3) ] <- stats::cov(t(scores[, , j]), t(scores[ , , k]),
                                                                           use = "pairwise.complete.obs")
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

# compute covariance matrix of ML estimation

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
#' @return A 4x4 matrix containing an estimate of
#' \eqn{ Cov( \theta_j, \theta_k)}.
#' @export
#'
#' @examples
#'
#' # generate data
#' cvrt <- (1:100)/100
#' d <- 16
#' coords <- matrix(20*abs(stats::rnorm(d*2)), ncol = 2)
#' x <- generateData(seed = 1, n = 100, temp.cov = cvrt, d= d, locations = coords)
#'
#' # compute ML estimate at each station
#' mle1 <- fit_scalegev(x[, 1], temp.cov =cvrt)
#' mle3 <- fit_scalegev(x[, 3], temp.cov =cvrt)
#'
#'
#' # compute covariance between ML estimate of station 1 and ML estimate of station 3
#' Gams <-  estimate_gammas(x[, c(1,3)], parmat = cbind(mle1$mle, mle3$mle),
#'                          temp.cov =cvrt, rel_par = FALSE)
#' compute_sigmajk(mle1$mle, mle3$mle, Gams[ 1:3, 4:6], cvrt,
#'                 solve(mle1$hessian), solve(mle3$hessian))
#'
compute_sigmajk <- function(par.j, par.k, Gammajk, temp.cov, Jinv.j , Jinv.k) {

  s11 <- Gammajk[1, 1]
  s12 <- Gammajk[1, 2]
  s13 <- Gammajk[1, 3]
  s21 <- Gammajk[2,1]
  s22 <- Gammajk[2,2]
  s23 <- Gammajk[2,3]
  s31 <- Gammajk[3,1]
  s32 <- Gammajk[3,2]
  s33 <- Gammajk[3,3]

  muj <- par.j[1]
  sigmaj <- par.j[2]
  gammaj <- par.j[3]
  alphaj <- par.j[4]

  muk <- par.k[1]
  sigmak <- par.k[2]
  gammak <- par.k[3]
  alphak <- par.k[4]

  fj <- (1-alphaj*temp.cov/muj)/sigmaj
  fk <- (1-alphak*temp.cov/muk)/sigmak

  gj <- alphaj*temp.cov/muj^2
  gk <- alphak*temp.cov/muk^2
  sig11 <- mean(fk * ( fj * s11 - gj*s12) -
                  gk*( fj*s12 - gj *s22), na.rm = TRUE) # alphaj*temp.cov*s22/muj^2), na.rm = TRUE)

  sig12 <-  mean((fj*s12 - gj*s22)/sigmak, na.rm = TRUE)

  sig13 <- mean(fj*s13 - gj*s23, na.rm = TRUE)

  sig14 <- mean( temp.cov*( ( fj*s11 - gj*s12 )/sigmak +  ( fj*s12 - gj*s22)/muk ), na.rm = TRUE)

  sig21 <- mean((s21*fk - s22*gk)/sigmaj, na.rm = TRUE)

  sig22 <- mean(s22/(sigmaj*sigmak), na.rm = TRUE)

  sig23 <- mean( s23/sigmaj, na.rm = TRUE)

  sig24 <- mean(temp.cov*(s21/(sigmaj*sigmak) + s22/(sigmaj*muk)), na.rm = TRUE)

  sig31 <- mean( s31*fk - gk*s32, na.rm = TRUE)

  sig32 <- mean(s32/sigmak, na.rm = TRUE)

  sig33 <- mean( s33, na.rm = TRUE)
  sig34 <- mean(temp.cov * ( s31/sigmak + s32/muk), na.rm = TRUE)

  sig41 <- mean(temp.cov *( (s11/sigmaj + s21/muj) *fk - gk * (s12/sigmaj + s22/muj) ), na.rm = TRUE)

  sig42 <- mean(temp.cov/sigmak*(s12/sigmaj + s22/muj), na.rm = TRUE)

  sig43 <- mean(temp.cov*(s13/sigmaj + s23/muj), na.rm = TRUE)

  sig44 <- mean(temp.cov^2*( (s11/sigmaj + s21/muj)/sigmak + (s12/sigmaj + s22/muj)/muk), na.rm = TRUE)

  Summands <- matrix( c(sig11, sig12, sig13, sig14,
                        sig21, sig22, sig23, sig24,
                        sig31, sig32, sig33, sig34,
                        sig41, sig42, sig43, sig44),
                      byrow = TRUE, nrow = 4)

  Jinv.j %*%  Summands %*% Jinv.k

}



compute_cov_pooled <- function(dat, par, temp.cov, hessmat) {


  n.dat <- nrow(dat)
  d <- ncol(dat)
  hessmatinv <- solve(d*hessmat) # multiply by d since hessian was computed on data transformed to vector
                                 # (i.e. it's an average of n*d observations instead of n)
  scores <- array(dim = c(3, n.dat, d))
  for( j in 1:d ) {
    scores[ , , j] <- score_standard_univ(dat[, j], theta = par, temp.cov = temp.cov, rel_par = FALSE)
  }

  score_sums <- apply(scores, c(1,2), sum)

  cov_score_sums <- stats::var(t(score_sums), na.rm = TRUE)

  Bs <- Bmat(par, temp.cov = temp.cov)
  Ts <- Tsigma_inv(par, temp.cov = temp.cov)

  summands <- purrr::map2(.x = Bs, .y = Ts, ~ {
    aa <- (.x %*% .y)
    return(aa %*% cov_score_sums %*% t(aa))
    }
    )

  (hessmatinv %*% ((Reduce("+", summands))/n.dat) %*% hessmatinv)/n.dat

}


compute_ci_rl <- function(rlhat, rl_var, level = .05) {
  quant_norm <- stats::qnorm(1-level/2)

  c("lower" = rlhat - quant_norm*sqrt(rl_var), "upper" = rlhat + quant_norm*sqrt(rl_var))
}
