compute_sigmajk_eff <- function(par.j, par.k, Gammajk, temp.cov, Jinv.j , Jinv.k) {

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

  fj <- (1-alphaj*temp.cov/muj)
  fk <- (1-alphak*temp.cov/muk)

  sig11 <- mean(fk/sigmak * ( fj * s11/sigmaj - alphaj*temp.cov/muj^2*s12) -
    alphak*temp.cov/muk^2*( fj*s12/sigmaj - alphaj*temp.cov*s22/muj^2), na.rm = TRUE)

  sig12 <-  mean((fj*s12/sigmaj - alphaj*temp.cov/muj^2*s22)/sigmak, na.rm = TRUE)

  sig13 <- mean(fj/sigmaj*s13 - alphaj*temp.cov/muj^2*s23, na.rm = TRUE)

  sig14 <- mean( temp.cov*( ( fj/sigmaj*s11 - alphaj*temp.cov/muj^2*s12 )/sigmak +  ( fj/sigmaj*s12 - alphaj*temp.cov/muj^2*s22)/muk ), na.rm = TRUE)

  sig21 <- mean((s21*fk/sigmak - s22*alphak*temp.cov/muk^2)/sigmaj, na.rm = TRUE)

  sig22 <- mean(s22/(sigmaj*sigmak), na.rm = TRUE)

  sig23 <- mean( s23/sigmaj, na.rm = TRUE)

  sig24 <- mean(temp.cov*(s21/(sigmaj*sigmak) + s22/(sigmaj*muk)), na.rm = TRUE)

  sig31 <- mean( s31*fk/sigmak - alphak*temp.cov/muk^2*s32, na.rm = TRUE)

  sig32 <- mean(s32/sigmak, na.rm = TRUE)

  sig33 <- mean( s33, na.rm = TRUE)
  sig34 <- mean(temp.cov * ( s31/sigmak + s32/muk), na.rm = TRUE)

  sig41 <- mean(temp.cov *( (s11/sigmaj + s21/muj) *fk/sigmak - alphak*temp.cov/muk^2 * (s12/sigmaj + s22/muj) ), na.rm = TRUE)

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
