## Prepare libraries
library(data.table)
library(dplyr)
library(R.utils)
library(progress)

## Prepare paths
path <- "D:/GBIF_INVERTEBRATE_RAW/DOWNLOAD_12-06-23_A-I/"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Identify terrestrial gastropods
# Read species list
sp_list <- fread(paste0(path, "Species_List/0022906-230530130749713.csv"))

# Get species
sp_list <- sp_list[grep("SPECIES", sp_list$taxonRank), ]
sp_list <- sp_list[, c("taxonKey","acceptedTaxonKey","acceptedScientificName",
                       "order","family","species","class")]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MOLLUSC DOWNLOAD ONLY ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get gastropod species
gastro_list <- sp_list[class == "Gastropoda", ]

# TaxonMatch tool only takes files of 1500 lines max
splits <- ceiling(nrow(gastro_list)/1500)
gsr_batch1 <- gastro_list[1:ceiling(nrow(gastro_list)/splits) , ]
gsr_batch2 <- gastro_list[(ceiling(nrow(gastro_list)/splits) + 1):(ceiling(nrow(gastro_list)/splits)*2) , ]
gsr_batch3 <- gastro_list[(ceiling(nrow(gastro_list)/splits)*2 + 1):nrow(gastro_list) , ]

# Save
#fwrite(gsr_batch1, paste0(path, "MolluscaBaseMatches/GBIF_species_as_input_b1.txt"), 
#       row.names = FALSE, col.names = FALSE, sep = "\t")
#fwrite(gsr_batch2, paste0(path, "MolluscaBaseMatches/GBIF_species_as_input_b2.txt"), 
#       row.names = FALSE, col.names = FALSE, sep = "\t")
#fwrite(gsr_batch3, paste0(path, "MolluscaBaseMatches/GBIF_species_as_input_b3.txt"), 
#       row.names = FALSE, col.names = FALSE, sep = "\t")

# Read MolluscaBase outputs
gsr_batch1 <- fread(paste0(path, "MolluscaBaseMatches/GBIF_species_as_input_b1_matched.txt"),
                    header = TRUE, sep = "\t", fill = TRUE)
gsr_batch2 <- fread(paste0(path, "MolluscaBaseMatches/GBIF_species_as_input_b2_matched.txt"),
                    header = TRUE, sep = "\t", fill = TRUE)
gsr_batch3 <- fread(paste0(path, "MolluscaBaseMatches/GBIF_species_as_input_b3_matched.txt"),
                    header = TRUE, sep = "\t", fill = TRUE)

# Re-format and identify unmatched records
gastro_list_matched <- rbind(gsr_batch1, gsr_batch2, gsr_batch3)
colnames(gastro_list_matched) <- c("taxonKey", "acceptedTaxonKey", "GBIF_acceptedScientificName",
                                   "Order", "Family", "Species", "Class", "AphiaID",
                                   "Match_Type", "MB_ScientificName", "AphiaID_accepted",
                                   "MB_acceptedScientificName", "isMarine",
                                   "isBrackish", "isFresh", "isTerrestrial")
table(gastro_list_matched$Match_Type)
# There are 46 unmatched species and 4 approximately matched (near_1, near_3, phonetic)

# We will check those manually using checklist bank for name matching
# And searching directly WORMS or MolluscaBase for the environment
to_check <- gastro_list_matched[Match_Type %in% c("", "near_1", "near_3", "phonetic"),
                                c("taxonKey", "acceptedTaxonKey", 
                                  "GBIF_acceptedScientificName", "MB_acceptedScientificName")]
#fwrite(to_check, paste0(path, "MolluscaBaseMatches/unmatched_taxa_MB_manual.csv"), 
#       row.names = FALSE)

# Re-read file after manual filling
checked_gastro <- fread(paste0(path, "MolluscaBaseMatches/unmatched_taxa_MB_manual.csv"), 
                         header = TRUE)

# For automatically matched taxa, a taxon is terrestrial if it has 1 in the isTerrestrial column
# For the manually matched taxa, a taxon is terrestrial if it has "terrestrial" as environment
terr_taxa_auto <- gastro_list_matched[isTerrestrial == 1, taxonKey]
terr_taxa_man <- checked_gastro[Environment_Manual == "terrestrial", taxonKey]
to_keep <- unique(c(terr_taxa_auto, terr_taxa_man))
to_remove <- gastro_list[!(taxonKey %in% to_keep), taxonKey]

# Filter species list and keep terrestrial gastropods only
sp_list <- sp_list[!(taxonKey %in% to_remove), ]

# Out of the terrestrial gastropods not matched between GBIF and MolluscaBase, one species seems to be
# a GBIF issue: Iberus gualtierianus is the same species as Iberus gueltieranus
# the GBIF name comes from NCBI where it is considered a synonym of I. gualteranus
# And BoL where I. gualteranus does not exist
# Moreover the verbatimScientificName for all occs of this taxonKey is I. gualterianus
# Which indicates that the name was wrongly filled by the provider
# We will correct this manually
correct_atk <- sp_list[species == "Iberus gualtieranus", acceptedTaxonKey]
correct_asn <- sp_list[species == "Iberus gualtieranus", acceptedScientificName]
sp_list[species == "Iberus gualtierianus", c("acceptedTaxonKey")] <- correct_atk
sp_list[species == "Iberus gualtierianus", c("acceptedScientificName")] <- correct_asn
sp_list[species == "Iberus gualtierianus", c("species")] <- "Iberus gualtieranus"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END OF GASTROPODA CLEANING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# And in the occurrence table
# Read GBIF records
all_occ <- fread(paste0(path, "Records/0022905-230530130749713.csv"))


# TO ADAPT TO GIVEN DOWNLOAD IF ISSUES THAT NEED CORRECTION
# Correct
correct_sn <- all_occ[which(all_occ$species == "Iberus gualtieranus"), "scientificName"][1]
all_occ[species == "Iberus gualtierianus", c("scientificName")] <- unlist(correct_sn)
all_occ[species == "Iberus gualtierianus", c("species")] <- "Iberus gualtieranus"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Split occurrence tables

# Remove duplicate columns
sp_list_match <- sp_list[ , -c("order", "family", "species", "class")]

# Merge records with species list
all_occ <- merge(all_occ, sp_list_match, by.x = "taxonKey", by.y="taxonKey")

# Number of accepted taxa
accepted_occ <- unique(all_occ$acceptedTaxonKey)

# The acceptedTaxon column can also refer to subspecies, we would like to work with species
subsp <- unique(sp_list$species)

# First we create the directories outside of the main loop

# Primary checks of the taxonomy tp avoid issues
nrow(sp_list[sp_list$class == "", ])
nrow(sp_list[sp_list$order == "", ])
nrow(sp_list[sp_list$family == "", ])
nrow(sp_list[sp_list$species == "", ])
# There are no issues with unknown class, order or family, however we will still account it in the script
# We will skip also any null taxa

if (nrow(sp_list[sp_list$species == "", ]) > 0){subsp <- subsp[-which(subsp == "")]}

# Create class directories
CLASSES <- unique(sp_list$class)
for (C in CLASSES){
  dir.create(paste0("D:/GBIF_INVERTEBRATE_RAW/GBIF_extracted_data/", C))
}

# Create order directories
ORDERS <- unique(sp_list$order)

for (O in ORDERS){
  CLASS <- sp_list[which(sp_list$order == O), ]
  CLASS <- as.character(CLASS[1, "class"])
  if (CLASS == ""){
    dir.create(paste0("D:/GBIF_INVERTEBRATE_RAW/GBIF_extracted_data/Unknown_class/", O))  
  }else{
    dir.create(paste0("D:/GBIF_INVERTEBRATE_RAW/GBIF_extracted_data/", CLASS,"/", O))
  }
}

# Create family directories
FAMILIES <- unique(sp_list$family)

for (f in FAMILIES){
  
  # if the family is unknown there may be many options for order or class, so we will deal this on the fly
  if(f == ""){
    next
    
  }else{
    ORDER <- sp_list[which(sp_list$family == f), ]
    ORDER <- as.character(ORDER[1, "order"])
    
    # We have the same issue with orders, the right class cannot be found
    if (ORDER == ""){
     next
    }
    
    CLASS <- sp_list[which(sp_list$order == ORDER), ]
    CLASS <- as.character(CLASS[1, "class"])
    if (CLASS == ""){
      dir.create(paste0("D:/GBIF_INVERTEBRATE_RAW/GBIF_extracted_data/Unknown_class/", ORDER, "/", f))
    } else {
      dir.create(paste0("D:/GBIF_INVERTEBRATE_RAW/GBIF_extracted_data/", CLASS, "/", ORDER, "/", f))
    }
  }
}

# Split the GBIF table
prbar <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed: :elapsedfull || Remaining: :eta]",
                          total = length(subsp),
                          complete = "=", incomplete = "-", current = ">", clear = FALSE, width = 100)

for (i in subsp){
  
  dirflag <- FALSE
  
  # Get occurrences
  occurrences <- all_occ[which(all_occ$species == i), ]
  
  # Get taxonomy
  FAMILY <- as.character(occurrences[1, "family"])
  ORDER <- as.character(occurrences[1, "order"])
  CLASS <- as.character(occurrences[1, "class"])
  
  # Account for taxonomy issues
  if (CLASS == ""){
    CLASS <- "Unknown_class"
  }
  if (ORDER == ""){
    ORDER <- "Unknown_order"
    dirflag <- TRUE
  }
  if (FAMILY == ""){
    FAMILY <- "Unknown_family"
    dirflag <- TRUE
  }
  
  # Create a new directory if needed
  if (dirflag){
    dir.create("D:/GBIF_INVERTEBRATE_RAW/GBIF_extracted_data/", CLASS, "/", ORDER, "/", FAMILY, "/",
               recursive = TRUE)
    cat("Created new directory: ", 
        "D:/GBIF_INVERTEBRATE_RAW/GBIF_extracted_data/", CLASS, "/", ORDER, "/", FAMILY, "/")
  }
  
  # Write the table
  fwrite(occurrences, 
         paste0("D:/GBIF_INVERTEBRATE_RAW/GBIF_extracted_data/",
                CLASS, "/", ORDER, "/", FAMILY, "/", i, ".csv"), 
         row.names = FALSE)
  
  prbar$tick()

}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Re-organize files appropriately
# Combine Bee and wasp downloads
# Split moth and butterfly downloads

# Bees
path_bees <- "D:/GBIF_INVERTEBRATE_RAW/GBIF_extracted_data/Insecta/Hymenoptera/"

# Create new folder
dir.create(paste0(path_bees, "Anthophila"))

# Move other bee family folders to new folder
bees <- c("Andrenidae", "Apidae", "Colletidae", "Halictidae", "Megachilidae", "Melittidae")
file.copy(from = paste0(path_bees, bees),
          to = paste0(path_bees, "Anthophila"),
          recursive = TRUE, copy.mode = TRUE, copy.date = TRUE)
file.remove(from = paste0(path_bees, bees))

# Same for wasps
path_wasps <- path_bees
dir.create(paste0(path_wasps, "Ichneumonoidea"))
wasps <- c("Braconidae", "Ichneumonidae")
file.copy(from = paste0(path_wasps, wasps),
          to = paste0(path_wasps, "Ichneumonoidea"),
          recursive = TRUE, copy.mode = TRUE, copy.date = TRUE)
file.remove(from = paste0(path_wasps, wasps))

# Butterflies
path_lepi <- "D:/GBIF_INVERTEBRATE_RAW/GBIF_extracted_data/Insecta/Lepidoptera/"
dir.create(paste0(path_lepi, "Rhopalocera"))
butterflies <- c("Hedylidae", "Hesperiidae", "Lycaenidae", "Nymphalidae",
                "Papilionidae", "Pieridae", "Riodinidae")
file.copy(from = paste0(path_lepi, butterflies),
          to = paste0(path_lepi, "Rhopalocera"),
          recursive = TRUE, copy.mode = TRUE, copy.date = TRUE)
file.remove(from = paste0(path_lepi, butterflies))

# Moths
dir.create(paste0(path_lepi, "Heterocera"))
lepi_families <- list.files(path_lepi, recursive = FALSE)
moths <- lepi_families[!lepi_families %in% c(butterflies, "Rhopalocera", "Heterocera")]
file.copy(from = paste0(path_lepi, moths),
          to = paste0(path_lepi, "Heterocera"),
          recursive = TRUE, copy.mode = TRUE, copy.date = TRUE)
file.remove(from = paste0(path_lepi, moths))
