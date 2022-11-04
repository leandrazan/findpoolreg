gev2frech <- function (x, loc, scale, shape) {
  if (shape == 0) {
    exp((x - loc)/scale)
  }
  else {
    pmax(1 + shape * (x - loc)/scale, 0)^(1/shape)
  }
}

frech2gev <- function (x, loc, scale, shape) {
  if (shape == 0) {
    scale * log(pmax(x, 0)) + loc
  }
  else {
    loc + scale * (pmax(x, 0)^shape - 1)/shape
  }
}


