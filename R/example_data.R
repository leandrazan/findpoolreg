#' Example data set
#'
#'
#' @format A 75 x 9 matrix of observations sampled from a scale-GEV model with spatial dependence structure
#' described by a max-stable process.
#'
#' @examples The data was generated with the following code.
#' meancoords <-  as.matrix(findpoolreg::example_grid[ , 6:7])
#' cvrt <- findpoolreg::GMST$smoothedGMST[68:142]
#' generateData(n = 75, d = 9,
#'             scale = c(rep(6, 3), rep(5, 3), rep(3, 3)),
#'             loc =  c(rep(24, 3), rep(20, 3), rep(17, 3)),
#'             shape =  c(rep(0.1, 6), rep(0.15, 3)),
#'             alpha =  c(rep(3, 3), rep(2.5, 3),  rep(2, 3)),
#'             seed = 1, locations = meancoords, temp.cov = cvrt)
#'
"example_data"





