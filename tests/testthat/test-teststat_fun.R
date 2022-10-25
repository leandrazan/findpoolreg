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


