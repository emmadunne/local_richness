countLocalRichness <- function(occs, rank = c("species","genus","indet","family"), per = "collection") {

	countIndetSpeciesRichness <- function(occs, taxonomic_hierarchies) {
		occs[grepl("informal", occs$identified_name), ]$accepted_name <- gsub(pattern = "informal", replacement = "", x = occs[grepl("informal", occs$identified_name), ]$identified_name)
		accepted_nos <- unique(occs$accepted_no) #need to search taxonomic hiearchy by accepted_number because some "accepted_names" now don't exist in the taxonomic hierarchy (e.g. "Chiroptera informal sp. 1")
		accepted_names <- unique(occs$accepted_name) #get all unique accepted taxonomic names
		X <- unique(unlist(taxonomic_hierarchies[match(accepted_nos, taxonomic_hierarchies$accepted_no), ][["t"]])) #get all the higher taxonomic levels for all occurrences
		sort(accepted_names[!(accepted_names %in% X)]) #get species plus indets (occurrences that don't represent higher taxonomic levels of other occurrences)
	}

	countIndetGenusRichness <- function(occs, taxonomic_hierarchies) {
		tmp <- occs[which(occs$accepted_rank %in% c("genus","species","subspecies") & !duplicated(occs$occurrence.genus_name)), c("accepted_name","occurrence.genus_name")] #get all occurrences assigned to genus or lower, removing duplicates
		definitelyuniquegenera_acceptednames <- unique(tmp$accepted_name) #get unique accepted names of unique genera (have to use full accepted names to search taxonomic hierarchy as some generic names don't have entries)
		definitelyuniquegenera <- unique(tmp$occurrence.genus_name) #get genus names that are definitely unique
		possiblyuniquegenera <- unique(occs[which(!(occs$accepted_rank %in% c("genus","species","subspecies"))), ]$accepted_name) #get names that are possibly unique genera
		X <- unique(unlist(taxonomic_hierarchies[match(c(possiblyuniquegenera,definitelyuniquegenera_acceptednames), taxonomic_hierarchies$accepted_name), ][["t"]])) #get combined taxonomic hierarchy (ranks above occurrence-name level)
		indet_genera <- possiblyuniquegenera[!(possiblyuniquegenera %in% X)] #get unique indet genera
		sort(unique(c(definitelyuniquegenera, indet_genera))) #combine with definitely unique genera
	}


	if (per == "collection") {

		coll_df <- tibble(collection_no = unique(occs$collection_no))

		if (rank == "species") {

			out <- distinct(occs[occs$accepted_rank %in% c("species","subspecies") & occs$occurrence.binomial != "", c('collection_no','occurrence.binomial')]) %>% group_by(collection_no) %>% summarise(species = list(sort(occurrence.binomial)))
			out <- left_join(coll_df, out)
			out <- out %>% mutate(species_richness = map(species, length)) %>% unnest(species_richness)

		} else if (rank == "genus") {

			coll_nos <- unique(occs$collection_no)
			genus_richness <- pbmclapply(coll_nos, function(x) {
				countIndetGenusRichness(occs = occs[which(occs$collection_no == x), ], taxonomic_hierarchies = taxonomic_hierarchies)
			}, mc.cores = detectCores(), mc.preschedule = TRUE)
			out <- data_frame(collection_no = coll_nos, genera = genus_richness)
			out <- out %>% mutate(genus_richness = map(genera, length)) %>% unnest(genus_richness)

		} else if (rank == "indet") {

			coll_nos <- unique(occs$collection_no)
			indet_richness <- pbmclapply(coll_nos, function(x) {
				countIndetSpeciesRichness(occs = occs[which(occs$collection_no == x), ], taxonomic_hierarchies = taxonomic_hierarchies)
			}, mc.cores = detectCores(), mc.preschedule = TRUE)
			out <- data_frame(collection_no = coll_nos, indets = indet_richness)
			out <- out %>% mutate(indet_richness = map(indets, length)) %>% unnest(indet_richness)

		} else if (rank == "family") {

			out <- tbl_df(table(distinct(occs[occs$accepted_rank %in% c("species","genus","subgenus","subspecies","family","subfamily") & occs$family != "", c('collection_no','family')])$collection_no))
			out <- distinct(occs[occs$accepted_rank %in% c("species","genus","subgenus","subspecies","family","subfamily") & occs$family != "", c('collection_no','family')]) %>% group_by(collection_no) %>% summarise(families = list(sort(family)))
			out <- left_join(coll_df, out)
			out <- out %>% mutate(family_richness = map(families, length)) %>% unnest(family_richness)

		}

	} else if (per == "formation") {

		fm_df <- tibble(formation = unique(occs$formation)); fm_df <- fm_df[fm_df != "", ]

		if (rank == "species") {

			out <- distinct(occs[occs$accepted_rank %in% c("species","subspecies") & occs$occurrence.binomial != "", c('formation','occurrence.binomial')]) %>% group_by(formation) %>% summarise(species = list(sort(occurrence.binomial)))
			out <- left_join(fm_df, out)
			out <- out %>% mutate(species_richness = map(species, length)) %>% unnest(species_richness)

		} else if (rank == "genus") {

			formations <- unique(occs$formation); formations <- formations[formations != ""]
			genus_richness <- pbmclapply(formations, function(x) {
				countIndetGenusRichness(occs = occs[which(occs$formation == x), ], taxonomic_hierarchies = taxonomic_hierarchies)
			}, mc.cores = detectCores(), mc.preschedule = TRUE)
			out <- data_frame(formation = formations, genera = genus_richness)
			out <- out %>% mutate(genus_richness = map(genera, length)) %>% unnest(genus_richness)

		} else if (rank == "indet") {

			formations <- unique(occs$formation); formations <- formations[formations != ""]
			indet_richness <- pbmclapply(formations, function(x) {
				countIndetSpeciesRichness(occs = occs[which(occs$formation == x), ], taxonomic_hierarchies = taxonomic_hierarchies)
			}, mc.cores = detectCores(), mc.preschedule = TRUE)
			out <- data_frame(formation = formations, indets = indet_richness)
			out <- out %>% mutate(indet_richness = map(indets, length)) %>% unnest(indet_richness)

		} else if (rank == "family") {

			out <- distinct(occs[occs$accepted_rank %in% c("species","genus","subgenus","subspecies","family","subfamily") & occs$family != "", c('formation','family')]) %>% group_by(formation) %>% summarise(families = list(sort(family)))
			out <- left_join(fm_df, out)
			out <- out %>% mutate(family_richness = map(families, length)) %>% unnest(family_richness)

		}

	}

	out[out[ ,3] == 0, ] <- NA
	return(out)
}
