# FUNCTIONS USED IN CLEANING SCRIPT P1

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Test if there are raws in a data frame
check_rows <- function(df){ifelse(nrow(df) > 0, TRUE, FALSE)}

# where df is a data.frame

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Sample randomly records on grid so that each grid cell only has 1 record and also retain record information
# Modified version of Plantarum's function from https://github.com/rspatial/dismo/issues/20
gridSampleTWS <- function (spdf, r, n = 1) {
  
  # Add record index
  spdf$record <- 1:length(spdf)
  
  # Extract Coordinates
  xy <- cbind(coordinates(spdf), index = spdf$record)
  
  # Extract cell number by coordinate
  r <- raster(r)
  cell <- cellFromXY(r, xy)
  
  # Warn if all points are outside the grid extent
  if (nrow(xy) == length(which(is.na(cell)))) {
    warning("All points outside grid extent")
    return(NA)}  
  
  # Numbers of non-NA cells
  uc <- unique(stats::na.omit(cell))
  
  # Add cell numbers to the table of coordinates 
  xy <- cbind(xy, cell = cell, rand = runif(nrow(xy)))
  
  # Remove missing cells
  xy <- stats::na.omit(xy)
  
  # seems unlikely there'll be unique rows, when one column
  # is a newly-generated random number?
  
  xy <- unique(xy)
  
  # Sort by random number to select observations randomly
  xy <- xy[order(xy[, "rand"]), ]
  
  # Create data frame and account for cases where one point is retained
  xy <- as.data.frame(xy)
  if (ncol(xy) == 1) {xy <- t(xy)}
  
  # Subset observations
  pts <- data.frame(numeric(), numeric(), numeric(),
                    numeric(), numeric())
  names(pts) <- names(xy)
  for (u in uc) {
    ss <- subset(xy, xy[, "cell"] == u)
    pts <- rbind(pts, ss[1:min(n, nrow(ss)), ])
  }
  
  # Ensure row names are preserved
  names(pts) <- colnames(xy)
  
  # Return result
  ret <- spdf[spdf$record %in% pts$index, ]
  ret <- ret[, names(ret) != "record"]
  return(ret)
}

# where spdf is a SpatialPointsDataFrame
# r is the raster defining the grid

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

standard_filters <- function(df, template_raster = NULL, 
                             coord_system = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs",
                             clean_coords = TRUE, clean_dataset = FALSE, on_grid = TRUE,
                             precision_col = "coordinateUncertaintyInMeters", precision = NULL,
                             lon = "decimalLongitude", lat = "decimalLatitude",
                             convert_codes = FALSE, countries = NULL,
                             outputs = c("table_prec", "summary", "final_rasters", "final_tables"), ...){
  
  # Final results
  result <- list()
  pre <- 0
  cc_coords <- 0
  cc_dataset <- 0
  grid_dupli <- 0
  end_obs <- 0
  
  # Start analysis
  start_obs <- nrow(df)
  
  # Filter precision
  if (!is.null(precision) & check_rows(df)){
    df <- df[df[ , precision_col] <= precision, ]
    pre <- start_obs - nrow(df)
  }
  
  # Convert country references to ISO3
  if(!is.null(countries) & convert_codes == TRUE & check_rows(df)){
    df[ , countries] <- countrycode(df[ , countries], origin =  'iso2c', destination = 'iso3c')
  }
  
  # Apply CoordinateCleaner coordinate tests
  if (clean_coords == TRUE & check_rows(df)){
    flags_cc <- suppressWarnings(clean_coordinates(x = df, lon = lon, lat = lat,
                                                   value = "spatialvalid",
                                                   verbose = FALSE, ...))
    
    cc_coords <- length(which(!flags_cc$.summary))
    df <- df[flags_cc$.summary, ]
  }
  
  # Apply CoordinateCleaner dataset tests
  if (clean_dataset == TRUE & check_rows(df)){
    df <- cd_round(cd_ddmm(df, lon = lon, lat = lat, value = "clean", verbose = FALSE, ...), 
                   lon = lon, lat = lat, value = "clean", verbose = FALSE, graphs = FALSE, ...)
    
    cc_dataset <- (start_obs - nrow(df)) - pre - cc_coords
  }
  
  # Save filtered observations
  if ("table_prec" %in% outputs){
    result[[length(result) + 1]] <- df
  }
  
  # Sample one observation by grid cell
  if (on_grid == TRUE & check_rows(df)){
    
    # Convert to spdf and check coordinate system
    df <- SpatialPointsDataFrame(df,
      coords = df[ , c(lon, lat)],
      proj4string = CRS(coord_system))
    
    if (coord_system != template_raster@crs@projargs){
      warning("Different CRS between observations and template: projecting...")
      df <- spTransform(df, template_raster@crs)
    }
    
    # Keep information of retained records with adapted grid sample
    if ("final_tables" %in% outputs){
      
      df <- gridSampleTWS(df, template_raster, n = 1)
      result[[length(result) + 1]] <- df
      if (suppressWarnings(is.na(df))){end_obs <- nrow(df@data)}
      grid_dupli <- (start_obs - end_obs) - (pre + cc_coords + cc_dataset)
      
    } else {
      
      # For faster and simpler we can also directly rasterize
      df <- raster::rasterize(df, template_raster, 1)
      end_obs <- cellStats(df, "sum")
      grid_dupli <- (start_obs - end_obs) - (pre + cc_coords + cc_dataset)
      
      if("final_rasters" %in% outputs) {
        
        result[[length(result) + 1]] <- df
        
      }
      
    }
    # Return an empty raster if no observations were retained during filtering
  } else if (on_grid == TRUE & !check_rows(df)){
    
    result[[length(result) + 1]] <- raster(template_raster)
    
  } else {
      
    end_obs <- nrow(df)
    
    }
  
  # Cleaning summary
  if ("summary" %in% outputs){
    
    result[[length(result) + 1]] <- c(start_obs, pre, cc_coords, cc_dataset, grid_dupli, end_obs)
    
  }
  
  return(result)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Calculate species richness during rasterization

spec_richness <- function(x) length(unique(x))

# where x is the vector of your data.frame/spdf containing the species names

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Combine species rasters

combine_rasters <- function(path_body, template, fun,
                            base = "D:/PostFiltering/"){
  
  # All possible data sources
  path_GBIF <- paste0(base, "GBIF_only/")
  path_OBSORG <- paste0(base, "OBSorg_only/")
  path_SG <- paste0(base, "SafeGuard/")
  path_GC <- paste0(base, "GlobalCollembola/")
  path_opCL <- paste0(base, "OpenEarthworms/")
  path_clCL <- paste0(base, "PrivateEarthworms/")
  
  paths <- c(paste0(path_GBIF, path_body, ".tif"),  # GBIF
             paste0(path_OBSORG, path_body, "_OBSORG.tif"), # Obsorg
             paste0(path_SG, path_body, "_SG.tif"), # SafeGuard
             paste0(path_GC, path_body, "_GC.tif"), # GlobalCollembola
             paste0(path_opCL, path_body, "_opCL.tif"), # Edaphobase + Phillips data
             paste0(path_clCL, path_body, "_clpCL.tif")) # Vers2022 + JM data
  
  # Empty raster to add observations
  species_raster <- template
  species_raster[] <- NA
  
  # Combine all rasters in 1
  for (path in paths){
    
    if (file.exists(path) & fun == "max"){
      species_raster <- max(species_raster, rast(path), na.rm = TRUE)
    } else if (file.exists(path) & fun == "sum"){
      species_raster <- sum(species_raster, rast(path), na.rm = TRUE)
    }
    
  }
  
  return(species_raster)
  
}