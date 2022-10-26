test_that("H0 LS functional works", {

  d <- sample(2:10, size = 1)
  A <- matrix(1: (4*d), ncol = d)

  expect_length(htheta(A, equal_distr = FALSE), 3*(d -1))

})

test_that("H0 ED functional works", {

  d <- sample(2:10, size = 1)
  A <- matrix(1: (4*d), ncol = d)

  expect_length(htheta(A, equal_distr = TRUE), 4*(d -1))

})

test_that("Derivative of H0 ED functional works", {

  d <- sample(2:10, size = 1)
  A <- matrix(1: (4*d), ncol = d)

  expect_equal(dim(grad_h(A, equal_distr = TRUE)) ,  c(4*(d-1), 4*d))
})


test_that("Derivative of H0 LS functional works", {

  d <- sample(2:10, size = 1)
  A <- matrix(1: (4*d), ncol = d)

  expect_equal(dim(grad_h(A, equal_distr = FALSE)) ,  c(3*(d-1), 4*d))
})


test_that("Singular covariance matrix produces warning", {
  xx <- matrix(rep(exp((1:100)/100), 3)*evd::rgev(300), ncol = 3)
  xx <- cbind(xx, xx[, 1])
  mlest <- fit_spat_scalegev(data = xx, temp.cov = (1:100)/100, hom = FALSE)
  expect_warning(teststat(mlest$mle, mlest$cov.mat, n = 100))
})


