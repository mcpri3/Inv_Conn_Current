
# Get the list of patrimonial invertebrates for France (from Patrinat)
myData <- openxlsx::read.xlsx(here::here('data/raw-data/SNAPTaxons_Representativite_Mlx_SpeciesLevel.xlsx'))
myData <- myData[myData$GROUP2_INPN %in% c('Arachnides', 'Bivalves', 'Crustacés', 'Entognathes', 'Gastéropodes', 'Insectes', 'Annélides', 'Myriapodes', 'Plathelminthes'),]
myData <- myData[!is.na(myData$LB_NOM_VALIDE_SPE_LEVEL), ]
myData <- dplyr::distinct(myData[, c('GROUP2_INPN', 'ORDRE', 'LB_NOM_VALIDE_SPE_LEVEL')])

table(myData$GROUP2_INPN)
by(myData$ORDRE, myData$GROUP2_INPN, table)

# GBIF check : Annélides, Plathelminthes ; no data 
# GBIF check : Myriapodes ; only 1 species with 23 occurrences over 49 species 
# GBIF check : Bivalves ; ok - 9 freshwater species 
# GBIF check : Arachnides ; only Araneae have occurrence data and 16/21 are doable
# GBIF check : Crustacés ; not much data + aquatic 
# GBIF check : Gastéropodes : Architaenioglossa and Stylommatophoraare terrestrial and have data available 
# GBIF check : Entognathes : Diplura no data and Collembola, to be check 
# GBIF check : Crustacés : no data / aquatic 

# Get the list of species in Marianne folder 
setwd("/Volumes/Seagate")
ffiles <- list.files('Invertebrates_CleanedOccur/GBIF_only/', recursive = T)
ffiles <- ffiles[grep('Species_rasters', ffiles)]   

sp.lst <- strsplit(ffiles, '/Species_rasters/')
sp.order <- lapply(sp.lst, function(x) {return(x[1])})
sp.lst <- lapply(sp.lst, function(x) {return(x[2])})
sp.lst <- lapply(sp.lst, function(x) {return(gsub('.tif', '', x))})
sp.lst <- data.frame(ORDER = unlist(sp.order), SPECIES = unlist(sp.lst))                 

# Compare the two 
idx <- myData$LB_NOM_VALIDE_SPE_LEVEL %in% sp.lst$SPECIES
in.it <- myData[idx, ]
by(myData$ORDRE, myData$GROUP2_INPN, table)
by(in.it$ORDRE, in.it$GROUP2_INPN, table)
out <- myData[!idx,]

# Check synonyms 
syn <- check_syn(out$LB_NOM_VALIDE_SPE_LEVEL, sp.lst$SPECIES)
syn$inDB <- syn$Synonym %in% sp.lst$SPECIES
syn <- syn[syn$inDB, ]

# Check occurrences
check.occur <- function(df, all.files = ffiles) {
  s <- df$LB_NOM_VALIDE_SPE_LEVEL
  f <- all.files[grep(s, all.files)]
  rr <- terra::rast(paste0('./Invertebrates_CleanedOccur/GBIF_only/', f))
  france <- sf::st_read(here::here('data/raw-data/SIG/Grids/CasestudyOutlines.gpkg'), layer = 'France')
  france <- sf::st_transform(france, sf::st_crs(rr))
  val <- terra::extract(rr, france)
  val <- sum(na.omit(val[, 2]))
  ddoable <- ifelse(val >= 5, 1, 0)
  return(data.frame(LB_NOM_VALIDE_SPE_LEVEL = s, doable = ddoable))
  }
library(dplyr)
check <- group_by(in.it, LB_NOM_VALIDE_SPE_LEVEL) %>% do(check.occur(df = .)) 
sum(check$doable)
myData <- left_join(myData, check, by = 'LB_NOM_VALIDE_SPE_LEVEL')
