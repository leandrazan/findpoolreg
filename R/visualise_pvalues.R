#' Visualise bootrstrap test results for pairs
#'
#' @param coord_grid A dataframe containing the grid coordinates in the columns `from_lon,to_lon, from_lat, to_lat`.
#' Can also contain the grid centers through columns `meanlon`, `meanlat`.
#' @param testres Result of bootstrap test as obtained from  \code{\link[findpoolreg]{bootstrap_scalegev_subsets}} or
#'  \code{\link[findpoolreg]{bootstrap_pairs_bivmod}}.
#' @param rej_proc The rejection procedure: either `holm` for the Holm stepdown procedure, `BY` for the
#' Benjamini Yekutieli stepup procedure, `BH` for the Benjamini Hochberg step-up procedure or `IM` when not adjusting for
#' multiple testing.
#' @param plot_type The plot type: either `rejection`, then regions are shaded according to test result,
#'  i.e. rejection or no rejection or `pvals`, then  regions are shaded according to their (adjusted) p-values.
#' @param level The test level
#' @param bins Bins for p-value shades when `plot_type = pvals`.
#' @param ... additional arguments passed to \code{\link[leaflet]{leafletOptions}}
#' @return A leaflet plot.
#' @export
#'
#' @examples
#'
#' ## generate data and compute test results
#' \dontrun{
#' coord_grid <- tidyr::expand_grid(from_lat = seq(44.75, 48.5, 1.25), from_lon = seq(7.25, 13.25, 2) ) %>%
#'          dplyr::mutate(to_lon = from_lon + 2, to_lat = from_lat + 1.25, Region = 1:16)
#'
#' meancoords <-   as.matrix((coord_grid %>% dplyr::mutate(meanlon = (from_lon + to_lon)/2,
#'                                                            meanlat = (from_lat + to_lat)/2))[ , 6:7])
#' cvrt <- (-20:79)/100
#' simdat <- generateData(n = 100, d = 16, scale = c(rep(6, 3), rep(5, 10), rep(3, 3)),
#'                        loc =  c(rep(24, 3), rep(20, 7), rep(17, 6)), shape =  c(rep(0.1, 13), rep(0.15, 3)),
#'                        alpha =  c(rep(3, 3), rep(2.5, 7),  rep(2, 6)),
#' seed = 1, locations = meancoords, temp.cov =  cvrt)
#'
#' subsets <- purrr::map(c(1:9, 11:16), ~ c(10, .x))
#' bootres <- bootstrap_scalegev_subsets(data = simdat, temp.cov = cvrt, locations = meancoords, subsets = subsets)
#'
#' plot_regions(coord_grid = coord_grid, testres = bootres, method = "holm", plot_type = "pvals")
#' plot_regions(coord_grid = coord_grid, testres = bootres, method = "BY")
#'
#' plot_regions(coord_grid = coord_grid, testres = bootres, plot_type = "pvals")
#'}
plot_regions <- function(coord_grid, testres, method = "holm", plot_type = "rejection", level = 0.1,
                         bins =  c(0,0.005, 0.01, 0.02, 0.05, 1), roi,  ...) {

  add.args <- list(...)

  d.roi <- length(roi)

  if(!("meanlon" %in% colnames(coord_grid))) {
   coord_grid <- coord_grid %>% dplyr::mutate(meanlon = (from_lon + to_lon)/2,
     meanlat = (from_lat + to_lat)/2)
  }

  testres <- get_adj_pvals(testres, methods = method, rejection = FALSE)

  if( "sbst" %in% colnames(testres)) {
    testres  <- testres %>% tidyr::unnest(cols = sbst)
  }
  testres <- testres %>% dplyr::rename("Region" = paste0("X", d.roi +1))


  testres <- testres %>% dplyr::rename( "padj" = dplyr::starts_with(paste0("p_", method), ignore.case = TRUE))

  testres <- testres %>% dplyr::mutate(fillop = ifelse(padj <= 0.1, 0.7, 0.5))

  if(plot_type == "rejection") {


   testres <- testres %>% dplyr::mutate(reject = (padj <= level))


   bb <- coord_grid %>% dplyr::left_join(testres, by = "Region")

   plotgrid <- leaflet::leaflet(data = bb, options = leaflet::leafletOptions( ... )) %>%
      leaflet::addTiles() %>%
      leaflet::addRectangles(lng1 = ~from_lon, lng2 = ~to_lon, lat1 = ~from_lat, lat2 = ~to_lat,
                    fill = ~reject, fillOpacity = 0.6) %>%
      leaflet::addCircles( lng = ~meanlon, lat = ~meanlat) %>%
      leaflet::addLabelOnlyMarkers(~meanlon, ~meanlat, label =  ~as.character(Region),
                          labelOptions =
                            leaflet::labelOptions(noHide = TRUE,
                                         textOnly = TRUE,
                                         textsize = 5,
                                         style = list(
                                           "color" = "blue",  "font-weight" = "bold",
                                           "font-size" = "30px")))
  }
  if(plot_type == "pvals") {

    bb <- coord_grid %>% dplyr::left_join(testres, by = "Region")

    mybins <- bins
    mypalette <-  leaflet::colorBin( palette = "YlGn", domain = (testres$p_boot), na.color="transparent",
                                     bins=mybins, reverse = TRUE)

    plotgrid <- leaflet::leaflet(data = bb) %>%
      leaflet::addTiles() %>%
      leaflet::addRectangles(lng1 = ~from_lon, lng2 = ~to_lon, lat1 = ~from_lat, lat2 = ~to_lat,
       fillColor  = ~mypalette(p_boot), fillOpacity =   ~fillop, color = "blue", fill = 1) %>%
      leaflet::addCircles( lng = ~meanlon, lat = ~meanlat) %>%
      leaflet::addLabelOnlyMarkers(~meanlon, ~meanlat, label =  ~as.character(Region),
                                   labelOptions =
                                     leaflet::labelOptions(noHide = TRUE,
                                                  textOnly = TRUE,
                                                  textsize = 5,
                                                  style = list(
                                                    "color" = "blue",  "font-weight" = "bold",
                                                    "font-size" = "30px"))) %>%
     leaflet::addLegend( pal=mypalette, values=~(p_boot),
                opacity=0.6, title = "(adjusted) p value",
                position = "bottomleft" )
  }

  plotgrid

}

