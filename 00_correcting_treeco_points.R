library(dplyr)
library(data.table)

##### ADJUSTING TreeCo DA/ha #####
points = readRDS('data/00_points_filtered.rds')
points = as.data.frame(points)
plots_coord <- points %>% dplyr::select(SiteCode.x, long1, lat1)

#Species list from abundance_data_set_per_sitecode_v1.csv
species_list = read.csv('treeco/treeco_for_guilherme_grittz/abundance_data_set_per_sitecode_v1.csv',
                        header = T, sep = ',')

#Filtering what is needed
#No filters for DA/ha, using everything — only filtering before IDW

#Removing some NA columns
species_list <- species_list[,-which(names(species_list) == "SubSite")]
species_list <- species_list[,-which(names(species_list) == "DoR")]
species_list <- species_list[,-which(names(species_list) == "family")]
species_list <- species_list[,-which(names(species_list) == "genus")]
sum(is.na(species_list))

#Final changes in species names
setDT(species_list)
species_list[species.correct %in% "Sebastiania brasiliensis", species.correct := "Sebastiania ramosissima"]
species_list[species.correct %in% "Lafoensia pacari", species.correct := "Lafoensia vandelliana"]
species_list[species.correct %in% "Prunus myrtifolia", species.correct := "Prunus brasiliensis"]
species_list[species.correct %in% "Pseudobombax grandiflorum", species.correct := "Pseudobombax majus"]
species_list[species.correct %in% "Randia armata", species.correct := "Randia ferox"]
species_list[species.correct %in% "Machaerium brasiliense", species.correct := "Machaerium paraguariense"]
species_list[species.correct %in% "Qualea cryptantha", species.correct := "Qualea glaziovii"]
species_list[species.correct %in% "Amaioua intermedia", species.correct := "Amaioua guianensis"]
species_list[species.correct %in% "Allophylus puberulus", species.correct := "Allophylus edulis"]
species_list[species.correct %in% "Piptocarpha organensis", species.correct := "Piptocarpha densifolia"]
species_list[species.correct %in% "Pleroma mutabile", species.correct := "Pleroma raddianum"]
species_list[species.correct %in% "Pterocarpus rohrii", species.correct := "Pterocarpus violaceus"]


#Correcting new synonyms from Flora do Brasil 2020 and stuff that occur in SC (IFFSC regional data)
syn.br <- read.csv("new_synonyms_floraBR.csv", na.strings = c(""," ",NA), as.is = TRUE)
syn.br <- syn.br[syn.br$status %in% c("replace", "invert"), ]
for (i in 1:dim(syn.br)[1]) {
  sp.i <- syn.br$original.search[i]
  rpl.i <- syn.br$search.str[i]
  st.i <- syn.br$status[i]
  
  if (st.i == "replace")
    species_list$species.correct[species_list$species.correct %in% sp.i] <- rpl.i
}

#Now IFFSC
syn.br <- read.csv("new_synonyms_IFFSC.csv", sep = ';')
syn.br <- syn.br[syn.br$status %in% c("replace", "invert"), ]
for (i in 1:dim(syn.br)[1]) {
  sp.i <- syn.br$original.search[i]
  rpl.i <- syn.br$search.str[i]
  st.i <- syn.br$status[i]
  
  if (st.i == "replace")
    species_list$species.correct[species_list$species.correct %in% sp.i] <- rpl.i
}

#Obtaining coordinates point for each register based on the SiteCode
joined_data <- merge(species_list, 
                     points, 
                     by.x = 'SiteCode', 
                     by.y = 'SiteCode.x', 
                     sort = FALSE, 
                     all.y = TRUE)

#Removing N less than 1
joined_data <- joined_data %>% filter(N.x >= 1)

#Re-measuring DA/ha for species (and converting all plots different from 1 ha sampling effort (e.g., 0.4 ha) to 1 ha)
#First, we re-measure N, since there are some differences in TreeCo
joined_data$Ntotal_new = with(joined_data, ave(N.x, SiteCode, FUN = sum))

#Then we adjust taking sampling effort into consideration
joined_data$DA_new = joined_data$Ntotal_new / joined_data$effort_ha.x

joined_data = joined_data[, c("SiteCode", "Ntotal_new", "DA_new")]
names(joined_data)[1] = "SiteCode.x"

points_newmetrics = merge(points, joined_data, by = "SiteCode.x", all.x = TRUE, sort = FALSE)
points_newmetrics = points_newmetrics %>%
  group_by(SiteCode.x) %>%
  slice(1)

saveRDS(points_newmetrics, 'data/00_points_filtered_nm.rds')
saveRDS(species_list, 'data/checked_specieslist.rds')

rm(list=ls())





