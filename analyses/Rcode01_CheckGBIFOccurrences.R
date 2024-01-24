#' @description 
#' 1. This script checks whether occurrence data are available (GBIF) for the pre-selected invertebrate groups and consequently narrows the list to non-data deficient 
#' species 
#' Filters applied : 
#' - Spatial coverage : metropolitan France 
#' - Temporal coverage : 2010 - 2023
#' - Minimal occurrence spatial precision : 1km 
#' - Type of occurrence : human observation
#' - Taxon level : species
#' - apply coordinate cleaner 
#' 
#' 2. This script generates distribution raster for each selected species 
#' 3. This script generates sampling effort raster for each selected order
#'  
 
#' @author Marie-Caroline Prima \email{marie-caroline.prima@univ-grenoble-alpes.fr}
#' 
#' @date 2024/01/08

# Get the list of patrimonial invertebrates for France (from Patrinat)
myData <- openxlsx::read.xlsx(here::here('data/raw-data/Taxa/SNAPTaxons_Representativite_Mlx_SpeciesLevel.xlsx'))
myData <- myData[myData$GROUP2_INPN %in% c('Insectes') | myData$ORDRE %in% c('Araneae', 'Architaenioglossa', 'Stylommatophora'),] #other groups are either aquatic or data deficient
myData <- myData[!is.na(myData$LB_NOM_VALIDE_SPE_LEVEL), ]
myData$LB_NOM_VALIDE_SPE_LEVEL[myData$LB_NOM_VALIDE_SPE_LEVEL == 'Prionotropis hystrix'] <- 'Prionotropis azami' #manual check of species syn. 
myData$LB_NOM_VALIDE_SPE_LEVEL[myData$LB_NOM_VALIDE_SPE_LEVEL == 'Omocestus navasi'] <- 'Omocestus antigai'
myData$LB_NOM_VALIDE_SPE_LEVEL[myData$LB_NOM_VALIDE_SPE_LEVEL == 'Bathysciola bonadonai'] <- 'Bathysciola brevicollis'
myData <- dplyr::distinct(myData[, c('GROUP2_INPN', 'ORDRE', 'CD_REF_SPE_LEVEL', 'LB_NOM_VALIDE_SPE_LEVEL')])
# Add family, class from Taxref 
taxref.full <- utils::read.csv(here::here('data/raw-data/Taxa/TAXREFv16.csv'), header = T, sep = ";")
taxref <- taxref.full[taxref.full$CD_REF %in% myData$CD_REF_SPE_LEVEL, ]
taxref <- dplyr::distinct(taxref[, c('CD_REF', 'PHYLUM', 'CLASSE', 'FAMILLE')])
myData <- dplyr::left_join(myData, taxref, by = c('CD_REF_SPE_LEVEL' = 'CD_REF'))

table(myData$GROUP2_INPN)
table(myData$ORDRE)
by(myData$ORDRE, myData$GROUP2_INPN, table)
by(myData$FAMILLE, myData$ORDRE, table)

# Get GBIF data, includes first round of cleaning 
all.files <- list.files(here::here('data/raw-data/GBIF/'), recursive = T)
all.files <- all.files[c(1, 2, 7, 8)]
all.occur <- data.table::data.table()
for (f in all.files) {
  ddta <- readr::read_delim(here::here(paste0('data/raw-data/GBIF/', f)),
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)
  ddta <- ddta[ddta$taxonRank == 'SPECIES',]
  ddta <- ddta[ddta$coordinateUncertaintyInMeters <= 1000,]
  ddta <- ddta[ddta$basisOfRecord %in% "HUMAN_OBSERVATION",]
  ddta <- ddta[ddta$year>= 2010, ]
  all.occur <- rbind(all.occur, ddta)
}

# Select on the same order that the species list
all.occur <- all.occur[all.occur$order %in% unique(myData$ORDRE),]
# Remove NA coordinates
all.occur <- all.occur[!is.na(all.occur$decimalLatitude), ]
## Clean occurrence data ##
all.occur <- CoordinateCleaner::cc_inst(all.occur) #remove presence in the vicinity of Biodiversity institutions
all.occur <- CoordinateCleaner::cc_val(all.occur)  #remove invalid coordinates 
all.occur <- CoordinateCleaner::cd_ddmm(all.occur,  ds = 'datasetKey', value = 'clean')  #remove presence from erroneous datasets 
# all.occur <- CoordinateCleaner::cd_round(all.occur, ds = 'datasetKey', graphs = F)  #remove presence from supposedly atlases 


#Read grid for spatial layers 
grid <- sf::st_read(here::here('data/raw-data/SIG/Grids/ReferenceGrid_France_bin_1000m.gpkg'))
grid.rast <- terra::rast(here::here('data/raw-data/SIG/Grids/ReferenceGrid_Europe_bin_1000m.tif'))
grid.rast <- terra::crop(grid.rast, grid)

# Project occurrences on the spatial grid and save spatial layers 
sp.df <- data.frame(species = myData$LB_NOM_VALIDE_SPE_LEVEL)
library(dplyr)
eval <- group_by(sp.df, species) %>% do(ProjectOccur(species = .)) %>% data.frame 

myData <- left_join(myData, eval, by = c('LB_NOM_VALIDE_SPE_LEVEL' = 'species'))
final <- myData[!is.na(myData$Ncell),]
final <- final[final$Ncell >= 5, ]
openxlsx::write.xlsx(final, './data/derived-data/PreFinalList-of-SNAP-Invertebrate-Species.xlsx')

# # Create the raster of sampling effort  
# for (g in unique(myData$ORDRE)) {
#   print(g)
#   sp.occur <- all.occur[all.occur$order %in% g,]
#   sp.occur <- sf::st_as_sf(sp.occur, coords = c("decimalLongitude" , "decimalLatitude"), crs = 4326) # Generate spatial layer of occurrences
#   sp.occur <- sf::st_transform(sp.occur, sf::st_crs(grid))
#   inter <- sf::st_intersects(grid, sp.occur)
#   inter <- lapply(inter, length)
#   grid$Nocc <- unlist(inter)
#   
#   sf::st_write(grid, here::here(paste0('data/derived-data/SIG/SamplingEffort/GBIFOccurrenceData_France_',
#                                         g,'_Res1000m_2010-2023.gpkg')), driver = 'GPKG', delete_layer = T)
#   
#   grid.rast <- terra::rasterize(grid, grid.rast, field = 'Nocc')
#   terra::writeRaster(grid.rast, here::here(paste0('data/derived-data/SIG/SamplingEffort/GBIFOccurrenceData_France_',
#                                                   g,'_Res1000m_2010-2023.tif')), overwrite = T) 
# 
# }

