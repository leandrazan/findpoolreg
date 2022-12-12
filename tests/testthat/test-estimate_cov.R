test_that("Efficient covariance estimation works", {
   cvrt <- (1:100)/100
   d <- 2
   coords <- matrix(20*abs(stats::rnorm(d*2)), ncol = 2)
   x <- generateData(seed = round(runif(1)*100), n = 100, temp.cov = cvrt, d= d, locations = coords)

   # compute ML estimate at each station
   par.j <- fit_scalegev(x[, 1], temp.cov =cvrt)
   par.k <- fit_scalegev(x[, 2], temp.cov =cvrt)

   Gams <-  estimate_gammas(x, parmat = cbind(par.j$mle, par.k$mle),
                            temp.cov =cvrt, rel_par = FALSE)
   Gammajk <- Gams[1:3, 4:6]

   s11 <- Gammajk[1, 1]
   s12 <- Gammajk[1, 2]
   s13 <- Gammajk[1, 3]
   s21 <- Gammajk[2,1]
   s22 <- Gammajk[2,2]
   s23 <- Gammajk[2,3]
   s31 <- Gammajk[3,1]
   s32 <- Gammajk[3,2]
   s33 <- Gammajk[3,3]

   muj <- par.j$mle[1]
   sigmaj <- par.j$mle[2]
   gammaj <- par.j$mle[3]
   alphaj <- par.j$mle[4]

   muk <- par.k$mle[1]
   sigmak <- par.k$mle[2]
   gammak <- par.k$mle[3]
   alphak <- par.k$mle[4]

   fj <- (1-alphaj*cvrt/muj)/sigmaj
   fk <- (1-alphak*cvrt/muk)/sigmak

   gj <- alphaj*cvrt/muj^2
   gk <- alphak*cvrt/muk^2

   sig11 <- fk * ( fj * s11 - gj*s12) -
                   gk*( fj*s12 - gj *s22)  # alphaj*temp.cov*s22/muj^2), na.rm = TRUE)

   sig12 <-  (fj*s12 - gj*s22)/sigmak

   sig13 <- fj*s13 - gj*s23

   sig14 <- cvrt*( ( fj*s11 - gj*s12 )/sigmak +  ( fj*s12 - gj*s22)/muk )

   sig21 <- (s21*fk - s22*gk)/sigmaj

   sig22 <-  s22/(sigmaj*sigmak)

   sig23 <-  s23/sigmaj

   sig24 <- cvrt*(s21/(sigmaj*sigmak) + s22/(sigmaj*muk))

   sig31 <-  s31*fk - gk*s32

   sig32 <- s32/sigmak

   sig33 <-  s33
   sig34 <- cvrt * ( s31/sigmak + s32/muk)

   sig41 <-cvrt *( (s11/sigmaj + s21/muj) *fk - gk * (s12/sigmaj + s22/muj) )

   sig42 <- cvrt/sigmak*(s12/sigmaj + s22/muj)

   sig43 <- cvrt*(s13/sigmaj + s23/muj)

   sig44 <- cvrt^2*( (s11/sigmaj + s21/muj)/sigmak + (s12/sigmaj + s22/muj)/muk)


   sigs <- list(sig11, sig21, sig31, sig41,
                sig12, sig22, sig32, sig42,
                sig13, sig23, sig33, sig43,
                sig14, sig24, sig34, sig44)

     Bj <- Bmat(par.j$mle, temp.cov = cvrt)
     Bk <- Bmat(par.k$mle, temp.cov = cvrt)

     Tj <- Tsigma_inv(par.j$mle, temp.cov = cvrt)
     Tk <- Tsigma_inv(par.k$mle, temp.cov = cvrt)

     BmalT.j <- purrr::map2(Bj, Tj, ~ {.x %*% .y})

     BmalT.k <- purrr::map2(Bk, Tk, ~ {.x %*% .y})

     Summands <- purrr::map2(BmalT.j, BmalT.k, ~ {.x %*% Gammajk %*% t(.y)} )

     for( j in 1:16) {

       expect_identical( round( purrr::map_dbl(Summands, ~ as.vector(.x)[j]) - sigs[[j]], 3), rep(0, 100) )

     }

})


test_that("Covariance estimation when NAs are present works", {
   cvrt <- (1:100)/100
   d <- 2
   coords <- matrix(20*abs(stats::rnorm(d*2)), ncol = 2)
   x <- generateData(seed = round(runif(1)*100), n = 100, temp.cov = cvrt, d= d, locations = coords)
   x[ c(3,6, 10), 1] <- NA

   # compute ML estimate at each station
   mle1 <- fit_scalegev(x[, 1], temp.cov =cvrt)
   mle2 <- fit_scalegev(x[, 2], temp.cov =cvrt)


   # compute covariance between ML estimate of station 1 and ML estimate of station 3
   Gams <-  estimate_gammas(x, parmat = cbind(mle1$mle, mle2$mle),
                            temp.cov =cvrt, rel_par = FALSE)

   testthat::expect_equal(dim(compute_sigmajk(mle1$mle, mle2$mle, Gams[ 1:3, 4:6], cvrt,
                   solve(mle1$hessian), solve(mle2$hessian))), c(4,4))
})
