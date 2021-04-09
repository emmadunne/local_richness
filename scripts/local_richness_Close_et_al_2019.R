# ******************************************************
#   
#   Local richness (alpha diversity)
#
# ______________________________________________________
#
#   Adapted by EM Dunne from:
#   Close et al. (2019) Nat. Ecol. Evol.
#
# ******************************************************


## Load packages (install first of course!)
library(tidyverse) # for data manipulation
library(pbmcapply) # for running various functions in parallel
library(stringr) # for running various functions
library(viridis) # cool colour schemes for plotting continuous data



# Getting set up ----------------------------------------------------------


## Data required (with example URLs for Late Triassic tetrapods):
##    A. PBDB occurrence data for your group and interval(s) of interest - e.g. http://paleobiodb.org/data1.2/occs/list.csv?datainfo&rowcount&base_name=Tetrapoda&ident=all&interval=Carnian,Rhaetian&private&show=full,classext,genus,ident,strat,env,ref,entname
##    B. PBDB collections data for the same group and interval(s) e.g. http://paleobiodb.org/data1.2/colls/list.csv?datainfo&rowcount&base_name=Tetrapoda&ident=all&interval=Carnian,Rhaetian&private&show=full,strat,lith
##    C. PBDB taxonomy data for the same group and interval(s) e.g. http://paleobiodb.org/data1.2/occs/taxa.csv?datainfo&rowcount&base_name=Tetrapoda&interval=Carnian,Rhaetian&private&show=full

## Files you'll need
##    i. This script (duh lol)
##    ii. The local richness algorithm file (countLocalRichness2019.R)

## Summary of steps:
##    1. Load and clean up the above datasets
##    2. Get the taxonomic hierarchies for the taxa in your occurrence dataset
##    3. Calculate local richness using Roger's algorithm (loaded from a separate R file)
##    4. Plot local richness using ggplot



# Load and clean data ---------------------------------------------------------------

## PBDB occurrence data
occ_data <- as_tibble(read.csv("./datasets/PBDB_occurrences_cleaned.csv", header = TRUE, stringsAsFactors=FALSE))

## Collection data download corresponding to the occurrence data above
coll_data <- as_tibble(read.csv("./datasets/PBDB_collections.csv", header = T, skip = 20, stringsAsFactors = FALSE))

## Taxonomic hierarchy information from the PBDB
tax_data <- as_tibble(read.csv("./datasets/PBDB_taxonomy.csv", stringsAsFactors = F, skip = 22, header = T))


## Make new fields: occurrence.genus_name and occurrence.species_name
occ_data[occ_data$accepted_rank == "species", "occurrence.genus_name"] <- word(occ_data[occ_data$accepted_rank == "species", ]$accepted_name, 1)
occ_data[occ_data$accepted_rank == "species", "occurrence.species_name"] <- word(occ_data[occ_data$accepted_rank == "species", ]$accepted_name, 2)
occ_data[occ_data$accepted_rank == "subspecies", "occurrence.genus_name"] <- word(occ_data[occ_data$accepted_rank == "subspecies", ]$accepted_name, 1)
occ_data[occ_data$accepted_rank == "subspecies", "occurrence.species_name"] <- word(occ_data[occ_data$accepted_rank == "subspecies", ]$accepted_name, 2)
occ_data[occ_data$accepted_rank == "genus", "occurrence.genus_name"] <- word(occ_data[occ_data$accepted_rank == "genus", ]$accepted_name, 1)
occ_data[occ_data$accepted_rank == "genus", "occurrence.species_name"] <- "sp."

## Make a binomial name field
occ_data$occurrence.binomial <- paste(occ_data$occurrence.genus_name, occ_data$occurrence.species_name, sep = " ")
occ_data[occ_data$occurrence.binomial == "NA NA", "occurrence.binomial"] <- ""

## Rename variables to match the ones Roger uses in his functions
occ_data$occurrence.reference_no <- occ_data$reference_no; occ_data$reference_no <- NULL
occ_data$occurrence.species_reso <- occ_data$species_reso; occ_data$species_reso <- NULL
occ_data$ma_max <- occ_data$max_ma; occ_data$max_ma <- NULL
occ_data$ma_min <- occ_data$min_ma; occ_data$min_ma <- NULL
occ_data$max_interval <- occ_data$early_interval; occ_data$early_interval <- NULL
occ_data$min_interval <- occ_data$late_interval; occ_data$late_interval <- NULL
occ_data$occurrence.abund_value <- occ_data$abund_value %>% as.numeric; occ_data$abund_value <- NULL
occ_data$country <- occ_data$cc
occ_data$paleolatdec <- occ_data$paleolat; occ_data$paleolngdec <- occ_data$paleolng

## Make other variables that will be useful for plotting
occ_data$ma_mid <- (occ_data$ma_max + occ_data$ma_min) / 2
occ_data$ma_length <- occ_data$ma_max - occ_data$ma_min

## Convert collection_no to a character in each dataset to avoid problems with merging
coll_data$collection_no <- as.character(coll_data$collection_no)
occ_data$collection_no <- as.character(occ_data$collection_no)

## Merge in useful variables from the collection data dataset
occ_data$collection.reference_no <- coll_data[match(occ_data$collection_no, coll_data$collection_no), ]$reference_no #slot in reference numbers for collections
occ_data$collection.ref_pubyr <- coll_data[match(occ_data$collection_no, coll_data$collection_no), ]$ref_pubyr #slot in ref_pubyrs for collections



# Taxonomic hierarchy -----------------------------------------------------

## This code grabs the taxonomic hierarchy for the species, genera, etc. in your occurrence dataset
## which is used later for calculating "indeterminate richness" (explained in next section)

## Function for taxonomic hierarchy
getTaxonomicHierarchy <- function(x) {
  x1 <- tax_data[which(tax_data$accepted_no == x), ]
  x2 <- vector()
  while (nrow(x1) != 0) {
    pn <- x1$parent_no
    x1 <- tax_data[which(tax_data$accepted_no == pn & tax_data$difference == ""), ]
    x2 <- unique(c(x1$accepted_name, x2))
  }
  x2
}

## Gather all the taxon numbers (accepted_no) for the taxa in the occurrence dataset
all_accepted_nos <- occ_data$accepted_no %>% unique

## Run the apply loop using mclapply 
system.time(taxonomic_hierarchies <- mclapply(all_accepted_nos, 
                                              function(x) getTaxonomicHierarchy(x), 
                                              mc.cores = detectCores(), 
                                              mc.preschedule = TRUE, 
                                              mc.cleanup = T)); names(taxonomic_hierarchies) <- all_accepted_nos

## Create a tibble for this output
taxonomic_hierarchies <- tibble(accepted_no = names(taxonomic_hierarchies), t = taxonomic_hierarchies)
taxonomic_hierarchies <- mutate(taxonomic_hierarchies, accepted_name = tax_data[match(taxonomic_hierarchies$accepted_no, tax_data$accepted_no), ]$accepted_name) %>% select(accepted_no, accepted_name, t)

## Some names may duplicate in the last step, this bit removes those so the richness-counting algorithm works properly
taxonomic_hierarchies <- pbmclapply(1:nrow(taxonomic_hierarchies), function(i) {
  taxa <- taxonomic_hierarchies[i, ]$t[[1]]
  taxonomic_hierarchies[i, ]$t[[1]] <- taxa[which(taxa != taxonomic_hierarchies[i, ]$accepted_name)]
  taxonomic_hierarchies[i, ]
  }, mc.cores = detectCores(), mc.cleanup = T, mc.preschedule = T)

## Clean up that output
taxonomic_hierarchies <- bind_rows(taxonomic_hierarchies)



# Calculating local richness ----------------------------------------------

### The main event - where we finally calculate local richness (or alpha diversity)

## Load the script containing richness-counting algorithm:
## *** be sure to change the file path to match where you've saved this file! ***
## (or open it as any normal script in R and click "Source" in the top right corner of the scripts pane)
source('~/Projects/_Sharing/local_richness/functions/countLocalRichness2019.R', echo=TRUE)


#### Calculate per-collection richness

#### So, species and genus richness are two measures that you don't necessarily need this method for
#### but it's the "indeterminate richness" or "IR" that we're really interested in 
#### (it's just neat to caculate the others alongside it in the same piece of code)
#### IR includes not only species richness, but also genera and other taxa that are indeterminate
#### at species-level, but that could feasibly represent a unique species
#### See methods section of Dunne et al. (2018) ProcB for an another description of this method
#### and Close et al. (2019) Nat Ecol Evol for the full method in action


## Start by truncating the dataset to just the columns the algorithm needs:
occs <- alpha_data %>% select(identified_name, difference, accepted_name, accepted_no, 
                              collection_no, accepted_rank, occurrence.binomial, 
                              occurrence.genus_name, family, formation)

### Now the calculations:

# Species richness
coll_SR <- countLocalRichness(occs, rank = "species", per = "collection")
colnames(coll_SR) <- c("collection_no", "genera", "SR") # rename the columns for ease later
# Genus richness
coll_GR <- countLocalRichness(occs, rank = "genus", per = "collection")
colnames(coll_GR) <- c("collection_no", "genera", "GR") # rename the columns for ease later
# Indet. richness
coll_IR <- countLocalRichness(occs, rank = "indet", per = "collection")
colnames(coll_IR) <- c("collection_no", "genera", "IR") # rename the columns for ease later


### Finally, organize the output and prepare for plotting

## Make a tibble for collection-level alpha diversity data from the coll_data imported at the beginning
collection_data <- coll_data[coll_data$collection_no %in% unique(alpha_data$collection_no), ]
collection_data <- left_join(collection_data, tbl_df(alpha_data[!duplicated(alpha_data$collection_no), c("collection_no","collection.reference_no","ma_max","ma_min","max_interval","min_interval","country","paleolatdec","paleolngdec","ma_mid","ma_length")]), by = "collection_no")

## Merge these data objects (richness data and collection data) into a single data frame for easier manipulation
collection_data <- left_join(collection_data, coll_SR, by = "collection_no")
collection_data <- left_join(collection_data, coll_GR, by = "collection_no")
collection_data <- left_join(collection_data, coll_IR, by = "collection_no")


## Subset and clean the output (collection_data) to use for plotting;
plotting_data <- subset(collection_data, select=c(collection_no, collection_name, cc, SR, GR, IR, paleolat, ma_mid, paleolatdec))
plotting_data <- na.omit(plotting_data) # remove any pesky NAs
plotting_data <- arrange(plotting_data, IR) # and finally, arrange this data from least to most rich in terms of IR




# Plotting in ggplot ------------------------------------------------------

## Define intervals - these are for the Late Triassic, so you'll need to change them to the ones you're plotting
interval_boundaries <- c(237, 228, 208.5, 201.3)

## Throw it all into ggplot:
alpha_plot <- ggplot(plotting_data, aes(x = ma_mid, y = paleolatdec, colour = IR)) +
  geom_vline(xintercept = interval_boundaries, lty = 2, col = "grey90") +
  geom_hline(yintercept = 0, colour = "grey10") +
  geom_point(alpha = 0.5, size = 7) +
  scale_color_viridis(trans = "log", breaks = c(1,5,22), direction = -1, option = "D") +
  scale_x_reverse(limits = c(237, 201), breaks = interval_boundaries) + 
  scale_y_continuous(labels = function(x) format(x, width = 5), limits = c(-90,90), breaks = seq(from = -90, to = 90, by = 15)) +
  theme_minimal() + theme(legend.position = c(0.2, 0.85), legend.direction = "vertical", panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), axis.title = element_text(size = 16)) +
  labs(x = "Time (Ma)", y = "Palaeolatitude (º)", colour = "Tetrapod local richness", cex = "Tetrapod local richness")
alpha_plot # open in RStudio plot window

## save the plot as a .pdf
ggsave("./plots/alpha_plot_IR.pdf", plot = alpha_plot, width = 28, height = 20, units = "cm")



