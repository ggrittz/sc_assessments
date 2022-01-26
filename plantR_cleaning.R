library(tidyverse)
library(plantR)
library(flora)

##################################################################
##### HERBARIA RECORDS FOR ALL TREE SPECIES THAT OCCUR IN SC #####
##################################################################


#Tree species that occur in SC (in prep. paper — names correct)
sc_trees_correctnames <- read.csv('sc_trees_list.csv', header = TRUE, sep = ';')
sc_trees_correctnames <- sc_trees_correctnames[order(sc_trees_correctnames$species), ]
sc_trees_correctnames <- sc_trees_correctnames$species


#GBIF points for SC tree species
plants_gbif <- rgbif2(species = sc_trees_correctnames,
                      country = "BR",
                      stateProvince = "Santa Catarina",
                      force = TRUE)
#plants_gbif <- readRDS('plants_gbif.rds')


#Multiple species from INCT needs a lapply
plants_inct <- lapply(sc_trees_correctnames, function(x) rspeciesLink(species = x,
                            basisOfRecord = "PreservedSpecimen",
                            Synonyms = "flora2020",
                            stateProvince = "Santa Catarina"))
names(plants_inct) <- sc_trees_correctnames
plants_inct <- bind_rows(plants_inct, .id = "original_search")

#plants_inct <- readRDS('plants_inct.rds')

#Preparing input data in the correct format
occs <- formatDwc(gbif_data = plants_gbif,
                  splink_data = plants_inct,
                  fix.encoding = c("splink_data", "gbif_data"),
                  drop = TRUE,
                  bind_data = TRUE)

#Format collection codes, people names, collector numbers, and dates
occs <- formatOcc(occs)

#Format locality information
occs <- formatLoc(occs)

#Format geographical coordinates
occs <- formatCoord(occs)

#Format species and family names
occs <- formatTax(occs)


##### DATA VALIDATION #####

#Validate localities
occs <- validateLoc(occs)

#Validate geographical coordinates
occs <- validateCoord(occs)
table(occs$geo.check)

#Validation species taxonomy and identification
occs <- validateTax(occs, generalist = TRUE, generalist.class = "medium") #add generalists


#Removing duplicates
occs$numTombo <- getTombo(occs[, "recordedBy.new"],
                          occs[, "recordNumber.new"])

#n_occur <- data.frame(table(occs$numTombo))
#n_occur[n_occur$Freq > 1,]
#occs[occs$numTombo %in% occs$numTombo[n_occur$Freq > 1],]

#Preparing to search for duplicates
dups <- prepDup(occs, comb.fields = list(c("family", "col.last.name", "col.number", "col.loc"),
                                         c("family", "col.last.name", "col.number", "col.year")))

dups <- getDup(dups)

occs <- cbind.data.frame(occs,
                         dups[, c("dup.ID", "dup.numb", "dup.prop")],
                         stringsAsFactors = FALSE)

occs <- mergeDup(occs)

table(occs$tax.check, occs$tax.check1)

occs <- rmDup(occs)

#Summary
summ <- summaryData(occs)

#Flags to take care
flags <- summaryFlags(occs)


##### FILTERING AFTER VALIDATION #####
#Filtering only records that occur in Brazil
occs1 <- occs %>% filter(country.new == "brazil")

#Filtering geographical coordinate
target1 <- c("ok_country", "ok_county", "ok_state", "ok_county_gazet", "ok_locality_gazet", "ok_state",
             "ok_state_gazet")
occs1 <- occs1 %>% filter(geo.check1 %in% target1)

#Filtering by determiner string
det_list <- read.csv('det_by_list.csv', header = TRUE, sep = ';')
occs1 <- occs1 %>% filter(identifiedBy.new1 %in% det_list$det_name)

#Removing cultivated species
occs1 <- occs1 %>% filter(is.na(cult.check))

#Changing subsp. and variety to species name only
occs1$suggestedName <- stringr::word(occs1$suggestedName, start = 1, end = 2)
sc_trees_correctnames$species <- stringr::word(sc_trees_correctnames$species, start = 1, end = 2)

#Removing any weird shit that may be left
occs1 <- occs1 %>% filter(suggestedName %in% sc_trees_correctnames)

write.csv(occs1, 'plantR_herbarium_filtered.csv')


