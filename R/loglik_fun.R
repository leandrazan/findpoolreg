
# params is a vector of length 4
# univariate fitting of scale model
# same as nll_univ_dj mit type = "scale" und rel_trend = FALSE

#' Negative log-likelihood of univariate GEV-scale-model
#'
#' @param params A vector of length 4
#' @param data The vector of observations
#' @param temp.cov Vector containing values of the temporal covariate.
#'
#' @return Numeric value of the log-likelihood of the data evaluated at the given parameters.
#' @export
#'
#' @examples
#' # generate some data
#' xx <- exp((1:100)/100)*evd::rgev(100)
#'
#' nll_scalegev(params = c(2,1, 0, 1), data = xx, temp.cov = (1:100)/100)
nll_scalegev <- function(params, data, temp.cov){

  if(!(length(data) == length(temp.cov))) { stop("'Data' and 'temp.cov' must have the same length.")}

  mu0 <- c("mu" = params[1])
  sigma0 <- c("sigma" = params[2])
  gamma0 <- c("gamma" = params[3])
  alpha0 <- c("alpha" = params[4])

  if(sigma0 <= 0) {return(1e+10)}
  else {
    mut <- mu0*exp(alpha0*temp.cov/mu0)
    sigmat <-  sigma0*exp(alpha0*temp.cov/mu0)


    if(abs(gamma0) < 1e-8){
      zt <- exp( -(data - mut)/sigmat)

      loglik <-  mean( log(sigmat) + log(zt) + zt , na.rm = TRUE)
    }
    else{
      zt <- 1 + gamma0*(data - mut)/sigmat
      if(any(zt < 0, na.rm = TRUE)){ loglik <- 1e+10}
      else{
        loglik <- mean(log(sigmat) + (1/gamma0 +1)*log(zt) + zt^(-1/gamma0), na.rm = TRUE)
      }
    }
    loglik
  }
}



#' Fit a scale-GEV-model to data
#'
#' @param data Either a numeric vector of observations, or a matrix containing
#'  observations from several stations/regions in its columns. In the latter case,
#' the prameters of the scale-GEV-model are assumed to be the same across columns,
#' such that only one set of parameters is computed.
#' @param temp.cov Values of the temporal covariate of the scale-GEV-model.
#' Must be of same record length as thr observations.
#' @param method Method passed to [stats::optim()].
#' @param maxiter Maximum number of iterations during maximisation (also passed to [stats::optim()]).
#' @param hessian logical; whether to return the numerically differentiated Hessian matrix.
#' @param printStartVals logical; whether to print start values of optimisation.
#'
#' @return A list containing the components
#' * mle: Estimated parameter values
#' * nll : Value of the negative log-likelihood evaluated at the optimised parameters
#' * conv : convergence code from [stats::optim()] (0 means successful convergence)
#' * hessian: The numerically differentiated Hessian matrix (when hessian was TRUE)
#' @export
#'
#' @examples
#' # generate univariate data
#' xx <- exp((1:100)/100)*evd::rgev(100)
#' fit_scalegev(data = xx, temp.cov = (1:100)/100, printStartVals = TRUE)
#'
#' # generate 3 dimensional data:
#' xx <- matrix(rep(exp((1:100)/100), 3)*evd::rgev(300), ncol = 3)
#' fit_scalegev(data = xx, temp.cov = (1:100)/100, printStartVals = TRUE)
fit_scalegev <- function(data, temp.cov, method = "BFGS", maxiter = 200,
                         hessian = TRUE, printStartVals = FALSE) {

  d.dat <- ncol(data)
  if(is.null(d.dat)) { d.dat <- 1}
  if(d.dat > 1) {
    n.dat <- nrow(data)
    # repeat values of temporal covariate for each data column
    if(length(temp.cov) == n.dat) { temp.cov <- rep(temp.cov, d.dat)}
    # transform data to one long vector
    data <- as.vector(data)
  }

  # compute start values from stationary model
  start_st <-  evd::fgev(as.vector(data), std.err = FALSE)$estimate
  start_vals <- c(start_st[1], start_st[2], start_st[3], 0)
  names(start_vals) <- c("mu", "sigma", "gamma", "alpha")
  if(printStartVals) {print(start_vals)}

  # optimise log-likelihood
  mlest <- stats::optim(start_vals, fn = nll_scalegev,
                 data = data, temp.cov = temp.cov,
                 method = method, control = list(maxit = maxiter),
                 hessian = hessian)

  if(!(mlest$convergence == 0) ){print("Optimization didn't succeed.")}

  if(!hessian) {
    return(list(mle = mlest$par, nll = mlest$value, conv = mlest$convergence))
  }
  else {

   list(mle = mlest$par, nll = mlest$value, conv = mlest$convergence,
                  hessian = mlest$hessian)
    }
}



# negative log Likelihood of spatial scale model under homogeneity constraint

#' Negative log-Likelihood of spatial scale-GEV-model under homogeneity constraint
#' @param params Vector of length \eqn{d +3}, where \eqn{d} is the number of sites.
#' Parameters must be given in the order
#' \eqn{(\mu_1, \mu_2, \ldots, \mu_d, \delta, \gamma, \eta)}.
#' @param data Matrix or data frame containing the observations from \eqn{d} measuring sites
#' (one column per site)
#' @param temp.cov Vector containing values of the temporal covariate.
#'
#' @return Numeric value of the log-Likelihood of the data evaluated at the given parameters.
#' @export
#'
#' @examples
#' # generate some data
#' xx <- (1:3)*matrix(rep(exp((1:100)/100), 3)*evd::rgev(300), ncol = 3)
#'
#' nll_scalegev_hom(params = c(1, 1, 1, 1, 0, 1), data = xx, temp.cov = (1:100)/100)
nll_scalegev_hom <- function(params, data, temp.cov){

  n.rec <- nrow(data)
  n.dat <- n.rec - colSums(is.na(data))
  d <- ncol(data)
  n.temp <- length(temp.cov)

  if(!(n.rec == n.temp)) { stop("Temporal covariate must have same length as record length.")}

  mu0s  <- params[1:d]
  names(mu0s) <- paste0("mu_", 1:d)
  delta <- c("delta" = params[d+1])

  gamma0 <- c("gamma" = params[d+2])
  # alpha0 <- c("alpha" = params[d+3])

  eta0 <- c("eta" = params[d+3])

  if(any(delta/mu0s <= 0)) {return(1e+10)}
  else {
    loglik <- numeric(d)
    for( j in 1:d) {
      mut <- mu0s[j]*exp(eta0*temp.cov)
      sigmat <-  mut/delta


      if(abs(gamma0) < 1e-8){
        zt <- exp( -(data[, j] - mut)/sigmat)

        loglik[j] <-  sum( log(sigmat) + log(zt) + zt , na.rm = TRUE)/n.dat[j]
      }
      else{
        zt <- 1 + gamma0*(data[,j] - mut)/sigmat
        if(any(zt < 0, na.rm = TRUE)){ loglik[j] <- 1e+10}
        else{
          loglik[j] <- sum(log(sigmat) + (1/gamma0 +1)*log(zt) + zt^(-1/gamma0), na.rm = TRUE)/n.dat[j]
        }
      }
    }

    sum(loglik)/d
  }
}


# sqrt(n) *(thatahat - theta)  \dto N(0, Sigma)
# sqrt(n) *(h(thetahat) - h(theta))  \dto N(0, h' Sigma h')
# unter H0: h(thetahat)  = 0
# sqrt(n)*h(theta)  \dto N(0, h' Sigma h')
# sqrt(n)*h(theta) (h' Sigma h')^(-1) sqrt(n)*h(theta) \approx chisq
# Sigmahat = cov(fisherinfoinv*score.fun)
# params is the vector (mu01, mu02, .., mu0d, delta, gamma, alpha)
# delta = mu0/sigma0 assumed constant across stations


# fit scale model to spatial data, possibly with homogeneity constraint.
# can return estimate of covariance matrix when not fitted with homogeneity assumption
#' Fit spatial scale-GEV-model to data
#'
#' @param data A matrix or dataframe representing the data, where each column corresponds to one site.
#' @param temp.cov Values of the temporal covariate of the scale-GEV-model.
#' Must be of same record length as thr observations.
#' @param hom logical; if `TRUE`, a spatial scale-GEV-model with homogeneity constraint is fitted.
#' If `FALSE`, a scale-GEV-model is fitted seperatly to each station. to fit model with homogeneity constraint. If
#' @param method  Method passed to [stats::optim()].
#' @param maxiter Maximum number of iterations during maximisation (also passed to [stats::optim()]).
#' @param returnRatios logical and only relevant if `hom = TRUE`: when TRUE, location-scale- and trend-location-parameter ratios are returned,
#' when FALSE, the plain parameters are returned.
#'
#' @return A list containing the estimated parameter values in `mle` and
#' * if `hom = TRUE`:  the value of the negative log-likelihood and a convergence code
#' * if `hom = FALSE`:  an estimation of the variance-covariance matrix.
#' @export
#'
#' @examples
#' # generate some data
#' xx <- (1:3)*matrix(rep(exp((1:100)/100), 3)*evd::rgev(300), ncol = 3)
#'
#' # fit componentwise scale-GEV-model
#' fit_spat_scalegev(data = xx, temp.cov = (1:100)/100, hom = FALSE)$mle
#'
#' # fit spatial scale-GEV-model with homogeneity constraint
#'
#' ## return global parameters, i.e. ratios of location/scale and trend/location
#' fit_spat_scalegev(data = xx, temp.cov = (1:100)/100, hom = TRUE, returnRatios = TRUE)$mle
#'
#' ## return local parameters, i.e. location, scale, shape and trend for each station
#' fit_spat_scalegev(data = xx, temp.cov = (1:100)/100, hom = TRUE, returnRatios = FALSE)$mle
#'
fit_spat_scalegev <- function(data, temp.cov, hom = FALSE, method = "BFGS",
                           maxiter = 200, returnRatios = TRUE) {

  if(!is.matrix(data)) { data <- as.matrix(data)}

  n.dat <- nrow(data)
  d <- ncol(data)


  if(hom) {

    start_st <-  fit_scalegev(data, temp.cov = temp.cov, hessian = FALSE)$mle
    start_vals <- c(rep(start_st[1], d), start_st[1]/start_st[2], start_st[3],
                    start_st[4]/start_st[1])

    names(start_vals ) <- c(paste0("mu_", 1:d), "delta", "gamma", "eta")
    # print(start_vals)
    mlest <- stats::optim(start_vals, fn = nll_scalegev_hom,
                   data = data, temp.cov = temp.cov,  method = method,
                   control = list(maxit = maxiter))
    if(!(mlest$convergence == 0) ){print("Optimization didn't succeed.")}
    mlpar <- mlest$par
    if(returnRatios) {
     mlpar <- matrix(c(mlpar[1:d], rep(mlpar[d+1],d),  rep(mlpar[d+2], d),
                      rep(mlpar[d+3], d)),
                    byrow = TRUE, nrow = 4)
     rownames(mlpar) <- c("mu", "delta", "gamma", "eta")
    } else {
      mlpar <- matrix(c(mlpar[1:d], mlpar[1:d]/rep(mlpar[d+1],d),  rep(mlpar[d+2], d),
                        rep(mlpar[d+3], d)*mlpar[1:d]),
                      byrow = TRUE, nrow = 4)
      rownames(mlpar) <- c("mu", "sigma", "gamma", "alpha")
    }
    return(list(mle = mlpar, nll = mlest$value, conv = mlest$convergence))
  }

  else {

    ml.est <- array(dim = c(4,d))

    Yfish <- array( dim = c(4*d, n.dat))
    for (m in 1:d) {

      estim <- fit_scalegev(data[, m], temp.cov = temp.cov, method = method,
                                maxiter = maxiter, hessian = TRUE)
      ml.est[, m] <- estim$mle
      score.m <- t(grad_ll_scalegev(data[, m], params = ml.est[,m],  temp.cov = temp.cov))
      #
      #       fishest <- cov(score.m, use = "pairwise.complete.obs")
      #       fishestinv <- solve(fishest)
      fishest <- estim$hessian
      fishestinv <- solve(fishest)  #, error = array(dim = c(4,4)))
      Yfish[(4 * m - 3):(4 * m), ] <- fishestinv %*% t(score.m)
    }

    cov.mat <- stats::var(t(Yfish), use = "pairwise.complete.obs")

    rownames(ml.est) <- c("mu", "sigma", "gamma", "aplha")




    return(list(mle = ml.est, cov.mat = cov.mat))
  }
}


#' Gradient of log-density in scale-GEV-model
#'
#' @param data Numeric vector of observations
#' @param params Vector of scale-GEV-model parameters (length 4)
#' @param temp.cov Values of temporal covariate (with same length as data)
#'
#' @return A \eqn{ 4\times} `length(data)`-matrix.
#' @export
#'
#' @examples
#' xx <- exp((1:100)/100)*evd::rgev(100)
#' grad_ll_scalegev(xx, params = c(1, 1, 0, 1), temp.cov = (1:100)/100)
grad_ll_scalegev <- function(data, params, temp.cov) {

  mu0 <- c(params[1])
  sigma0 <- c(params[2])
  xi <- c(params[3])
  alpha0 <- c(params[4])

  mut <- mu0*exp(alpha0*temp.cov/mu0)
  sigmat <-  sigma0*exp(alpha0*temp.cov/mu0)


  zt <- (data - mut)/sigmat
  in_supp <- which(1+ xi*zt > 0)
  zt[!in_supp] <- NA

  if(abs(xi) < 1e-08) {
    ut <- exp(-zt)
  }
  else {
    ut <- (1+ xi *zt)^(-1/xi)
  }

  dmu0ut <- -ut^(xi +1)*(alpha0 *temp.cov *data/(mu0^2*sigmat) - 1/sigma0)

  # ut^(xi +1)/sigmat* (expo*(1 - alpha0*temp.cov/mu0) - (x - mut)*alpha0*temp.cov/mu0)
  dsigma0ut <- - ut^(xi +1)/sigma0*(-data/sigmat + mu0/sigma0)

  dalpha0ut <- ut^(xi +1)/sigmat* data*temp.cov/mu0

  scoreloc0 <- alpha0*temp.cov/mu0^2 + dmu0ut*( (xi +1)/ut -1)

  scorescale <- -1/sigma0 + dsigma0ut*( (xi +1)/ut -1)

  scorealpha <- -temp.cov/mu0  + dalpha0ut*( (xi +1)/ut  -1)

  if(abs(xi) < 1e-08){
    scoreshape <- (1-ut)*zt^2/2- zt
  }
  else{
    scoreshape <- (1-ut)*1/xi*(1/xi*log(1+xi*zt)-zt/(1+xi*zt)) - zt/(1+xi*zt)
  }

  return(matrix(c(scoreloc0, scorescale, scoreshape, scorealpha),
                nrow = 4, byrow = TRUE))
}




