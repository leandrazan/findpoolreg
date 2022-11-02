test_that("Generating unit Frechet bootstrap samples works", {

  d <- sample(5:10, size = 1)
  coords <-  matrix(20*abs(stats::rnorm(d*2)), ncol = 2)
  xx <- generateData(seed = sample(1:1000, 1),
                     d = d, locations = coords,
                     scale = c(1:d), loc = rep(20, d), shape = seq(0, 0.3, length.out = d), alpha = rep(2.5, d))

  bootsamps <- generate_bootsamp_unitfrech(xx, temp.cov = (1:100)/100, locations = coords,
                                           ms_models = c("powexp", "gauss", "brown"), B = 4, sel_crit = "AIC", warn_msfit = FALSE)
  expect_identical(dim(bootsamps$X_star[[1]]), c(100L, d))
  expect_length(bootsamps$X_star, 4)
})
