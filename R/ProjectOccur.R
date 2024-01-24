#' ProjectOccur : project the number of occurrences per pixels on a provided gridded geopackage and raster 
#'
#' @param species single species name, obtained from the dply::group_by function applied to a data frame containing species list identified by 'species' column 
#' @param occur data frame of occurrences including X, Y coordinates and species columns  
#' @param colnams.coords X and Y column name for coordinates in occur data frame
#' @param coords.ref.sys reference system for provided coordinates
#' @param grid.gpkg sf object of the grid
#' @param grid.raster raster of the grid
#' @param nlim minimum number of occupied pixels for the spatial layers to be saved 
#'
#' @return Return a data frame with species evaluation : Total number of occurrences (Nocc) and number of occupied pixels (Ncell). Export the grid (if nlim is respected) with the number of 
#' occurrences in each pixel in geopackage and raster format (into folder /data/derived-data/SIG/Occurrence)
#' @export 
#'
#' @examples
ProjectOccur <- function(species, occur = all.occur, colnams.coords =  c("decimalLongitude", "decimalLatitude"), coords.ref.sys = 4326, grid.gpkg = grid, 
                         grid.raster = grid.rast, nlim = 5) {
  s <- species$species
  sp.occur <- occur[occur$species %in% s, ]
  Nocc <- nrow(sp.occur)
  
  if (Nocc >= nlim) {
    
    sp.occur <- sf::st_as_sf(sp.occur, coords = colnams.coords, crs = coords.ref.sys) # EPSG:4326 / Lat-Long
    sp.occur <- sf::st_transform(sp.occur, sf::st_crs(grid.gpkg))
    inter <- sf::st_intersects(grid.gpkg, sp.occur)
    inter <- lapply(inter, length)
    grid.gpkg$Nocc <- unlist(inter)
    Ncell <- sum(grid.gpkg$Nocc!=0) 
    
    if (sum(grid.gpkg$Nocc!=0) >= nlim) {
      
      # Save geopackage
      sf::st_write(grid.gpkg, here::here(paste0('data/derived-data/SIG/Occurrence/GBIFOccurrenceData_France_',
                                                gsub(' ','_', s),'_Res1000m_2010-2023.gpkg')), driver = 'GPKG', delete_layer = T)
      
      # Save raster
      grid.raster <- terra::rasterize(grid.gpkg, grid.raster, field = 'Nocc')
      terra::writeRaster(grid.raster, here::here(paste0('data/derived-data/SIG/Occurrence/GBIFOccurrenceData_France_',
                                                        gsub(' ','_',s),'_Res1000m_2010-2023.tif')), overwrite = T)
    } 
  } else {
    Ncell <- NA
  }
  return(data.frame(species = s, Nocc = Nocc, Ncell = Ncell))
}
