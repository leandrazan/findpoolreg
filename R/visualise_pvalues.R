## to do: different step up procedures, different adjusted p values
## visualising p value
## compare resutls of boottrap_scalegev_subsets and bootstrap_scalegev

#
# # visualise p vals and rejections
# bb <- tidyr::expand_grid(from_lat = seq(44.75, 53.5, 1.25), from_lon = seq(-0.75, 13.25, 2) ) %>%
#   dplyr::mutate(to_lon = from_lon + 2, to_lat = from_lat + 1.25, Region = 1:64, meanlon = (from_lon + to_lon)/2,
#                 meanlat = (from_lat + to_lat)/2)

plot_regions <- function(coord_grid, testres) {

  coord_grid <- coord_grid %>% dplyr::mutate(meanlon = (from_lon + to_lon)/2,
                           meanlat = (from_lat + to_lat)/2)

  testres  <- testres %>% tidyr::unnest(cols = sbst)

  bb <- coord_grid %>% dplyr::left_join(testres)

  leaflet::leaflet(data = bb) %>%

    leaflet::addTiles() %>%
    # addCircles(~from_lon, ~from_lat, radius=900, stroke=T) %>%
    leaflet::addRectangles(lng1 = ~from_lon, lng2 = ~to_lon, lat1 = ~from_lat, lat2 = ~to_lat,
                  fill = ~reject, fillOpacity = 0.6) %>%
    # color  = ~mypalette(p_corr), fillOpacity =  0.5 ) %>%
    # fillOpacity =  ~ (1 - p_corr), fillColor = "green") %>%
    leaflet::addCircles( lng = ~meanlon, lat = ~meanlat,
    ) %>%
    #  addRectangles(lng1 =  5.250  ,  lng2 = 7.250, lat1= 48.500 ,  lat2 = 49.750 )%>%
    #  addRectangles(lng1 = 7.250 ,  lng2 =  9.250, lat1= 49.750 , lat2=  51.000) %>%
    leaflet::addLabelOnlyMarkers(~meanlon, ~meanlat, label =  ~as.character(Region),
                        labelOptions =
                          labelOptions(noHide = TRUE,
                                       textOnly = TRUE,
                                       textsize = 5,
                                       style = list(
                                         "color" = "blue",  "font-weight" = "bold",
                                         "font-size" = "20px")))

}

#
# bb %>% left_join(rejyn)
#
# mybins <-  c(0,0.005, 0.01, 0.02, 0.05,  1)
# mypalette <- colorBin( palette = "YlGn", domain = (rejres2er$p_corr), na.color="transparent", bins=mybins,
#                        reverse = TRUE)
#
# leaflet(data = bb %>% left_join(rejyn)) %>%
#
#   addTiles() %>%
#   # addCircles(~from_lon, ~from_lat, radius=900, stroke=T) %>%
#   addRectangles(lng1 = ~from_lon, lng2 = ~to_lon, lat1 = ~from_lat, lat2 = ~to_lat,
#                 fill = ~reject, fillOpacity = 0.6) %>%
#   # color  = ~mypalette(p_corr), fillOpacity =  0.5 ) %>%
#   # fillOpacity =  ~ (1 - p_corr), fillColor = "green") %>%
#   addCircles( lng = ~meanlon, lat = ~meanlat,
#   ) %>%
#   #  addRectangles(lng1 =  5.250  ,  lng2 = 7.250, lat1= 48.500 ,  lat2 = 49.750 )%>%
#   #  addRectangles(lng1 = 7.250 ,  lng2 =  9.250, lat1= 49.750 , lat2=  51.000) %>%
#   addLabelOnlyMarkers(~meanlon, ~meanlat, label =  ~as.character(Region),
#                       labelOptions =
#                         labelOptions(noHide = TRUE,
#                                      textOnly = TRUE,
#                                      textsize = 5,
#                                      style = list(
#                                        "color" = "blue",  "font-weight" = "bold",
#                                        "font-size" = "20px")))
#
#
# ### schnecke on
