test_that("Transformation of scale-GEV to unit Frechet works", {

  xx <- evd::rfrechet(100)
  yy <- 5*exp(2/10*(1:100)/100)*(xx-1) + 10*exp(2/10*(1:100)/100)

  xxtest <- scalegev2frech(data = yy, temp.cov = (1:100)/100, par = c(10, 5, 1, 2))

  expect_equal(xx, xxtest)

})

test_that("d-dimensional transformation of scale-GEV to unit Frechet works", {

  xx <- generateData(seed = sample(1:1000, 1),
                     d = 3,
                     scale = c(1,2, 3), loc = rep(20, 3), shape = c(0, 0.2, 0.3, 9), alpha = rep(2.5, 9))

  yy <- spat_scalegev2frech(xx, temp.cov = (1:100)/100)

  for( i in 1:ncol(xx)) {
    xx[, i] <- scalegev2frech(xx[ ,i], temp.cov = (1:100)/100)
  }

  expect_identical(xx, yy)

})

