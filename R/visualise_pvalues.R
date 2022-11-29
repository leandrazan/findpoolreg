#' Visualise bootstrap test results for pairwise testing
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
#' @param loi Location(s) of interest, passed as numeric vector containing location labels. Can be omitted when location of interest is of length 1.
#' @param ... additional arguments passed to \code{\link[leaflet]{leafletOptions}}, e.g. `zoomControl = FALSE`. Further,
#' the arguments `position` (one of `bottomright, bottomleft, topright, topleft`) which gives the legend position, `width` and `height`
#' for the plot size can be given.
#' @return A leaflet plot.
#' @export
#'
#' @examples
#'
#' \dontrun{
#' ## generate data and compute test results
#' coord_grid <- tidyr::expand_grid(from_lat = seq(44.75, 47.25, 1.25), from_lon = seq(7.25, 11.25, 2) ) %>%
#'          dplyr::mutate(to_lon = from_lon + 2, to_lat = from_lat + 1.25, Region = 1:9)
#'
#' meancoords <-   as.matrix((coord_grid %>% dplyr::mutate(meanlon = (from_lon + to_lon)/2,
#'                                                            meanlat = (from_lat + to_lat)/2))[ , 6:7])
#' cvrt <- (-20:79)/100
#' simdat <- generateData(n = 100, d = 9, scale = c(rep(6, 3), rep(5, 3), rep(3, 3)),
#'                        loc =  c(rep(24, 3), rep(20, 3), rep(17, 3)), shape =  c(rep(0.1, 6), rep(0.15, 3)),
#'                        alpha =  c(rep(3, 3), rep(2.5, 3),  rep(2, 3)),
#' seed = 1, locations = meancoords, temp.cov =  cvrt)
#'
#' subsets <- purrr::map(c(1:4, 6:9), ~ c(5, .x))
#' bootres <- bootstrap_scalegev_subsets(data = simdat, temp.cov = cvrt, locations = meancoords, subsets = subsets, B = 100)
#'
#' visualise_test_res(coord_grid = coord_grid, testres = bootres, method = "holm", plot_type = "pvals", loi = 5)
#' visualise_test_res(coord_grid = coord_grid, testres = bootres, method = "BY", loi = 5)
#'
#' pvisualise_test_res(coord_grid = coord_grid, testres = bootres, method = "BH", plot_type = "pvals", loi = 5)
#'}
visualise_test_res <- function(coord_grid, testres, method = "holm", plot_type = "rejection", level = 0.1,
                         bins =  c(0, 0.05, 0.075,  0.1, 1), loi = NULL,  ...) {

  add.args <- list(...)
  if(is.null(add.args$position)) {add.args$position <- "bottomleft"}
  if(is.null(add.args$width)) {add.args$width <- 800}
  if(is.null(add.args$height)) {add.args$height <- 400}
  d.loi <- length(loi)
  if(is.null(loi)) { d.loi <- 1 }

  if(!("meanlon" %in% colnames(coord_grid))) {
   coord_grid <- coord_grid %>% dplyr::mutate(meanlon = (from_lon + to_lon)/2,
     meanlat = (from_lat + to_lat)/2)
  }

  testres <- get_adj_pvals(testres, methods = method, rejection = FALSE)

  if( "sbst" %in% colnames(testres)) {
    testres  <- testres %>% tidyr::unnest(cols = sbst)
  }
  testres <- testres %>% dplyr::rename("Region" = paste0("X", d.loi + 1))

  if(is.character(testres$Region)) {
    coord_grid$Region <- as.character(coord_grid$Region)
  }

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
    mypalette <- leaflet::colorBin(palette = "Reds",
                                   domain = seq(0,1, 0.05), na.color = "transparent",
                                   bins = mybins, reverse = TRUE)

    plotgrid <- leaflet::leaflet(data = bb, options = leaflet::leafletOptions( ... ),
                                 width = add.args$width, height = add.args$height) %>%
      leaflet::addTiles() %>%
      leaflet::addRectangles(lng1 = ~from_lon, lng2 = ~to_lon, lat1 = ~from_lat, lat2 = ~to_lat,
       fillColor  = ~mypalette(padj), fillOpacity =   ~fillop, color = "blue", opacity  = 1) %>%
      leaflet::addCircles( lng = ~meanlon, lat = ~meanlat, color = "blue") %>%
      leaflet::addLabelOnlyMarkers(~meanlon, ~meanlat, label =  ~as.character(Region),
                                   labelOptions =
                                     leaflet::labelOptions(noHide = TRUE,
                                                  textOnly = TRUE,
                                                  textsize = 5,
                                                  style = list(
                                                    "color" = "blue",  "font-weight" = "bold",
                                                    "font-size" = "30px"))) %>%
     leaflet::addLegend( pal=mypalette, values=~(padj),
                opacity=0.6, title = "(adjusted) p value",
                position =  add.args$position )
  }

  plotgrid

}

