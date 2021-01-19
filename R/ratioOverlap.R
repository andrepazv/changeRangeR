#' @title Calculate the proportion of a range area that is either 1: contained by landcover categories, or 2: correlated with
#' a continuous environmental layer.
#' @description Calculate the proportion of the species' range (e.g., a thresholded SDM) that is contained by landcover categories
#' taken from a shapefile. Example shapefile categories include protected areas, threatened areas. ratioOverlap returns an s4 object
#' containing the masked raster layer and the percent of the total range that lies within the shapefile polygons specified.
#' @param r either raster or shapefile object representing a binary range.
#' @param shp (optional) either 1) a shapefile of land cover features or 2) a continuousnraster. Must be in same projection as r parameter. If shp is a raster, then the number of cells within each quantile are calculated.
#' @param rasMask (optional) a raster layer to calculate the Pearson correlation with the object r. Only if r or shp is a raster layer
#' @param field The shapefile field attribute containing the features to compare (i.e., the column name).
#' @param category a list of the names of shapefile features to compare. If all features are to be used, input "All".
#' @param proj character string proj4string of crs of landcover layer.
#' @export
#'
#'



ratioOverlap <- function(r, shp = NULL, rasMask = NULL, field, category){
  #setClass("ratioOverlap", slots = list(maskedRange = "RasterLayer", ratio = "character"))
  require(sf)
  require(rgdal)
  require(raster)
  require(rgeos)
  require(dplyr)

## if r is a shapefile
if(class(r) != "RasterLayer" & class(shp) != "RasterLayer"){
  r <- rgeos::gBuffer(r, byid = T, width = 0)
  shp <- rgeos::gBuffer(shp, byid = T, width = 0)
  r <- sf::st_as_sf(r)
  shp <- sf::st_as_sf(shp)

  if(category == "All"){
   maskedRange <- sf::st_intersection(r, shp)
  }else{
    fc <- lapply(category, function(x) dplyr::filter(shp, shp[[field]]==x))
    fc <- do.call("rbind", fc)
    maskedRange <- sf::st_intersection(r, fc)
  }
  ratio <- (sf::st_area(maskedRange) / sf::st_area(r)) * 100
  ratio <- paste0("Percentage of range within shape is ", ratio, "%")
  correlation <- NULL
}else{
  if(class(r) != "RasterLayer" & class(shp) == "RasterLayer"){
    r <- raster::rasterize(r, shp, method = "ngb")
    maskedRange <- raster::mask(r, shp)
    maskedQuants <- raster::quantile(raster::mask(shp, r))
    q25 <- raster::ncell(shp[shp < maskedQuants[2]])
    q50 <- raster::ncell(shp[shp > maskedQuants[2] & shp < maskedQuants[3]])
    q75 <- raster::ncell(shp[shp > maskedQuants[3] & shp < maskedQuants[4]])
    ratq25 <- q25 / raster::ncell(r[!is.na(shp)])
    ratq50 <- q50 / raster::ncell(r[!is.na(shp)])
    ratq100 <- q100 / raster::ncell(r[!is.na(shp)])
    ratio <- rbind(paste0("The proportion of the range below 25%: ", ratq25), paste0("The proportion of the range between 25% and 50%: ", ratq50), paste0("The proportion of the range between 50% and 75%: ", ratq75),
                   paste0("The proportion of the range between 75% and 100%: ", ratq100))

    if(!is.null(rasMask)){
      rasMask.resam <- raster::resample(rasMask, r, method = "bilinear")
      rRasmask <- raster::stack(rasMask.resam, r)
      layerCorrs <- raster::layerStats(rRasmask, stat = "pearson")
      correlation <- layerCorrs$`pearson correlation coefficient`[[2]]
    }
  }
}


## if r is a raster
if(class(r)=="RasterLayer"){
  if(class(shp)[1] == "SpatialPolygonsDataFrame" | class(shp)[1] == "sf"){
    if (category == "All"){
      shp <- sf::st_as_sf(shp)
      r <- raster::projectRaster(r, crs = crs(shp)@projargs, method = 'ngb')
      maskedRange <- raster::mask(r, shp)
    }
    else {
      shp <- sf::st_as_sf(shp)
      # if (crs == NULL){
      #   fc <- filter(shp, grepl(category, field))
      #   out<- mask(r, fc)
      # } else {
      r <- raster::projectRaster(r, crs = crs(shp))
      fc <- lapply(category, function(x) dplyr::filter(shp, shp[[field]]==x))
      fc <- do.call("rbind", fc)
      #fc <- filter(shp, shp[[field]]==category)
      maskedRange <- raster::mask(r, fc)
      #  }
      #  return(out)
    }
    ratio <- raster::ncell(maskedRange[!is.na(maskedRange)]) / raster::ncell(r[!is.na(r)]) * 100
    ratio <- paste0("Percentage of range within shape is ", ratio, "%")
  }
  if(class(shp) == "RasterLayer"){
    shp <- raster::crop(shp, r)
    r <- raster::crop(r, shp)
    maskedRange <- raster::mask(r, shp)
    maskedQuants <- raster::quantile(raster::mask(shp, r))
    q25 <- raster::ncell(shp[shp < maskedQuants[2]])
    q50 <- raster::ncell(shp[shp > maskedQuants[2] & shp < maskedQuants[3]])
    q75 <- raster::ncell(shp[shp > maskedQuants[3] & shp < maskedQuants[4]])
    q100 <- raster::ncell(shp[shp > maskedQuants[4] & shp < maskedQuants[5]])
    ratq25 <- q25 / raster::ncell(r[!is.na(shp)])
    ratq50 <- q50 / raster::ncell(r[!is.na(shp)])
    ratq75 <- q75 / raster::ncell(r[!is.na(shp)])
    ratq100 <- q100 / raster::ncell(r[!is.na(shp)])
    ratio <- rbind(paste0("The proportion of the range below 25%: ", ratq25), paste0("The proportion of the range between 25% and 50%: ", ratq50), paste0("The proportion of the range between 50% and 75%: ", ratq75),
                   paste0("The proportion of the range between 75% and 100%: ", ratq100))
  }
  correlation = NULL
  if(!is.null(rasMask)){
    rasMask.resam <- raster::resample(rasMask, r, method = "bilinear")
    rRasmask <- raster::stack(rasMask.resam, r)
    layerCorrs <- raster::layerStats(rRasmask, stat = "pearson")
    correlation <- layerCorrs$`pearson correlation coefficient`[[2]]
  }
}

  #out <- new("ratioOverlap", maskedRange = maskedRange, ratio = ratio)
  return(list(maskedRange = maskedRange, ratio = ratio, correlation = correlation))

}

