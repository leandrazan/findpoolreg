
# computes \dot \ell_{(0, 1, \gamma)}(( x- mut)/sigmat)
score_standard_univ <- function(x, theta, temp.cov = NULL, rel_par = TRUE) {

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


estimate_gammas <- function(data, parmat, temp.cov = NULL, rel_par = TRUE) {
  n <- nrow(data)
  d <- ncol(data)

  scores <- array(dim = c(3, n, d))
  Gammas <- array(dim = c(3*d,3*d))
  for( j in 1:d ) {
    scores[ , , j] <- score_standard_univ(data[, j], theta = parmat[, j], temp.cov = temp.cov, rel_par = rel_par)
  }
  for (j in 1:d) {
    for (k in j:d) {
      Gammas[ ((j -1)*3 +1) : (j*3), ((k -1)*3 +1) : (k*3) ] <-
        Gammas[ ((k -1)*3 +1) : (k*3), ((j -1)*3 +1) : (j*3) ] <- cov(t(scores[, , j]), t(scores[ , , k]))
    }
  }
  Gammas
}

Bmat <- function(par, temp.cov) {

  n <- length(temp.cov)

  mu <- par[1]
  sigma <- par[2]
  alpha <- par[4]

  expo <- exp(alpha/mu * temp.cov)

  b11 <- ( 1- alpha*temp.cov/mu)*temp.cov
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

Tsigma_inv <- function(par, temp.cov) {

  mu <- par[1]
  sigma <- par[2]
  alpha <- par[4]
  sigmat <- sigma*exp(alpha/mu*temp.cov)

  purrr::map(sigmat, ~ diag(c(1/.x, 1/.x, 1) ))

}


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

