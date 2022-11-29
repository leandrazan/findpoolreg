#' Compute RL of GEV distribution in reference climate
#'
#' @param theta Estimated GEV parameters
#' @param Tyrl The return periods for which the RLs are computed
#' @param type The model: 'stationary'. 'shift' or 'scale'
#' @param ref_gmst The reference value of the temporal covariate for which the
#' RL is computed. Ignored when \code{type = "stationary"}.
#' @return A tibble containg the estimated RL for each combination of \code{Tyrl} and
#' \code{ref_gmst}.
#' @export
#'
#' @examples
#' compute_rl(theta = c(mu = 20, sigma = 2, gamma = 0.2, alpha = 2),
#' Tyrl = c(50, 100), type = "scale", ref_gmst = c(0.5, 0.95))
#'
compute_rl <- function(theta, Tyrl, type, ref_gmst = NULL) {


  if(type == "stationary") {
    rl <- evd::qgev(1-1/Tyrl, loc = theta[1], scale = theta[2], shape = theta[3])
    return(dplyr::tibble(rl = rl, Year = Tyrl ))
  }
  if(type == "scale") {

      mut <- theta[1]*exp(theta[4]/theta[1]*ref_gmst)
      sigmat <- theta[2]*exp(theta[4]/theta[1]*ref_gmst)

    return(purrr::map_dfr(Tyrl, ~ {
      rl <- evd::qgev(1-1/.x, loc = mut, scale = sigmat, shape = theta[3])
      dplyr::tibble(rl = rl, refGMST = ref_gmst, Year = .x)
    }))
  }
  if(type == "shift") {
    mut <- theta[1] + theta[4]*ref_gmst
    return(purrr::map_dfr(Tyrl, ~ {
      rl <- evd::qgev(1-1/.x, loc = mut, scale = theta[2], shape = theta[3])
      dplyr::tibble(rl = rl, refGMST = ref_gmst, Year = .x)
    }))
  }

}
