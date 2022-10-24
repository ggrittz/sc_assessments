rm(list=ls())

library(rredlist)
library(dplyr)
library(flora)

#List of species in SC
spp = readRDS('data/threat_mean_pop_sizes_for_ConR.rds')

#colClean <- function(x){colnames(x) <- gsub("trees_", "", colnames(x)); x} 
#spp = colClean(spp)

#Remover NAs
spp = na.omit(spp)

spp = spp[, 1]
spp = as.data.frame(spp)

##### Getting the IUCN assessments for Brazil (CNCFlora) - National level #####
tmp = flora::get.taxa(spp$spp, replace.synonyms = FALSE, life.form = TRUE)
tmp$search.str = spp$spp
tmp <- merge(spp, tmp, by.x = "spp", by.y = "search.str", all.x = TRUE, sort = FALSE)
#table(spp$spp == tmp$spp)
spp$status.reflora <- tmp$threat.status

##### Getting the global IUCN assessments (IUCN) ####
#Citation: IUCN 2020. IUCN Red List of Threatened Species. Version 2020-1 <www.iucnredlist.org>
iucn <- read.csv("IUCN_2020_assessments.csv", as.is = TRUE, na.strin = c(NA,""," "))

iucn$species.correct2 <- sapply(strsplit(iucn$scientificName," "), function(x) paste(x[1], x[2],sep=" "))

tmp <- merge(spp, iucn[, c("species.correct2" ,"redlistCategory", "redlistCriteria")], 
             by.x = "spp", by.y = "species.correct2", all.x = TRUE)

tmp <- tmp[,c("spp", "status.reflora", "redlistCategory", "redlistCriteria")]
table(spp$spp == tmp$spp)

##### GET LOCAL ASSESSMENT FROM CONSEMA/2014 #####
consema <- read.csv("consema_2014.csv", as.is = TRUE, header = TRUE, sep = ';')
consema_spp <- consema$species
consema_check <- flora::get.taxa(consema_spp, replace.synonyms = TRUE, life.form = TRUE)
consema$species <- consema_check$search.str
table(consema$species == consema_check$search.str)

#Adding species not found in CNCFlora (CONSEMA mistakes)
consema[270, 2] <- "Dendrophorbium missionum"
consema[42, 2] <- "Heleocharis montevidensis"
consema[98, 2] <- "Apoclada diversa"
consema[129, 2] <- "Elaphoglossum jamesonii"

final <- merge(tmp, consema[, c("species", "consema")], 
               by.x = "spp", by.y = "species", all.x = TRUE, sort = FALSE)
final <- final[, c(1, 5, 2, 3, 4)]

saveRDS(final, "rds_files/previous_assessments_spp.rds")

rm(list=ls())



