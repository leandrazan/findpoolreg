## ---- include = FALSE-----------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup----------------------------------------------------------------------------------------------------------------------------------
library(findpoolreg)
library(dplyr)

## ---- echo = FALSE--------------------------------------------------------------------------------------------------------------------------
data.frame(Location = 1:9, mu = c(23, 23, 22, 20, 20, 20, 20, 20, 21), 
           sigma = c(6, 6, 6, 5, 5, 6, 5, 5, 6), 
           gamma = c( 0.12, 0.1, 0.15, 0.1, 0.1, 0.15, 0.11, 0.1, 0.05), 
           alpha = c(rep(3, 3), rep(2.5, 6)))

## -------------------------------------------------------------------------------------------------------------------------------------------
data("example_grid")
data("example_data")
head(example_grid)
head(example_data)

## -------------------------------------------------------------------------------------------------------------------------------------------
leaflet::leaflet(data = example_grid, width = 500, height = 400)  %>% 
  leaflet::addTiles() %>%
  leaflet::addRectangles(lng1 = ~from_lon, lng2 = ~to_lon, lat1 = ~from_lat, lat2 = ~to_lat) %>% 
  leaflet::addLabelOnlyMarkers(~meanlon, ~meanlat, label =  ~as.character(Region), 
                               labelOptions = leaflet::labelOptions(noHide = TRUE, textOnly = TRUE,
                               textsize = 5, style = list("color" = "blue",  "font-weight" = "bold", 
                                                          "font-size" = "30px")))

## -------------------------------------------------------------------------------------------------------------------------------------------
colnames(example_data)
subsets <- purrr::map(c(1:4, 6:9), ~ c(5, .x))
subsets[1:2]  # first two subsets
cvrt <- (findpoolreg::GMST %>% dplyr::filter(Year >= 1947, Year <= 2021))$smoothedGMST

## -------------------------------------------------------------------------------------------------------------------------------------------
head(example_grid)
grid_centers <- as.matrix(example_grid[ , c("meanlon", "meanlat")])

