
test_that("gradient is correct", {
  xx <- evd::rgev(1, shape = 0.2)
  a <- -c(numDeriv::grad( function(x) { nll_scalegev(params = c(x, 2, 0.3, 0.4 ),
                                           data = xx,  temp.cov = 0.5) },
                  x = 1),
          numDeriv::grad( function(x) { nll_scalegev(params = c(1, x, 0.3, 0.4 ),
                                                     data = xx,  temp.cov = 0.5) },
                          x = 2),
          numDeriv::grad( function(x) { nll_scalegev(params = c(1, 2, x, 0.4 ),
                                                     data = xx,  temp.cov = 0.5) },
                          x = 0.3),
          numDeriv::grad( function(x) { nll_scalegev(params = c(1, 2, 0.3, x ),
                                                     data = xx,  temp.cov = 0.5) },
                                                  x = 0.4)
  )
  b <- as.vector(grad_ll_scalegev(data = xx, params = c(1,2,0.2, 0.4), temp.cov = 0.5))
  expect_equal(round(a-b), rep(0,4))
  })


