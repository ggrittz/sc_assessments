##### COMBINING ASSESSMENTS FROM ALL CRITERIA - A, B, C, D #####

rm(list=ls())

library(ConR)
library(data.table)
library(tidyverse)

#####################
#### PRE-EDITING ####
#####################
#Loading and merging the assessments
#Criterion A
critA <- readRDS("data/criterionA_all_GLs.rds")
critA <- critA[, c("species", "assessment.period", "reduction_A12", "A1","A2", "category_A",
                  "category_A_code", "reduction_A12.25ys", "A1.25ys", "A2.25ys")]
critA$reduction_A12 <- round(as.double(critA$reduction_A12),2)
critA$reduction_A12.25ys <- round(as.double(critA$reduction_A12.25ys),2)
  
#Criterion B
critB_opt <- readRDS("data/criterionB_herb.rds")
critB_opt <- critB_opt[, c("tax", "EOO", "AOO", "Nbe_subPop", "nbe_loc_total", "protected", "declineB",
                          "sever.frag", "category_B1", "category_B2", "category_B", "category_B_code")]

#Criterion C
critC <- readRDS("data/criterionC_all_prop_mature.rds")
critC <- critC[, c("species", "any.decline", "cont.decline", "C1", "C2", 
                   "category_C", "category_C_code", "C1.p0.45")]

#Criterion D
critD <- readRDS("data/criterionD_all_prop_mature.rds")
critD <- critD[, c("species", "pop.size", "pop.size.low", "D", "D.AOO", 
                   "D2.Loc", "category_D", "category_D_code", "D.p0.45")]
critD$pop.size <- round(as.double(critD$pop.size), 1)
critD$pop.size.low <- round(as.double(critD$pop.size.low), 1)

#Estimated pop. sizes based on AOO, taxonomy, life-form and endemism level
#est.pop <- readRDS("data/estimated_pop_size_from_AOO.rds")
#critD <- merge(critD, est.pop[, c("species", "pred", "pred.low")], by = "species", all = TRUE, sort = FALSE)
#critD$pop.size[is.na(critD$pop.size)] <- round(as.double(critD$pred[is.na(critD$pop.size)]), 1)
critD$pop.size.low[is.na(critD$pop.size.low)] <- round(as.double(critD$pred.low[is.na(critD$pop.size.low)]), 1)
critD <- critD[order(critD$species), ]
critD <- critD[, c("species", "pop.size", "pop.size.low", "D", "D.AOO",
                  "D2.Loc", "category_D", "category_D_code")]


##### Merging all assessments #####
all.crit <- merge(critA, critB_opt, by.x = "species", by.y = "tax", all = TRUE)
all.crit <- merge(all.crit, critC, by = "species", all = TRUE)
all.crit <- merge(all.crit, critD, by = "species", all = TRUE)
all.crit <- all.crit[order(all.crit$species), ]
rm(critA, critB_opt, critC, critD)

#Remove notata
#all.crit <- all.crit[-c(453), ]
all.crit <- all.crit %>% dplyr::rename(B1 = category_B1, B2 = category_B2)

subcriteria <- c("A1", "A2", "B1", "B2", "C1", "C2", "D")
for(i in 1:length(subcriteria)) {
  all.crit[, subcriteria[i]] <- 
    ConR::near.threatened(cats = all.crit[ , subcriteria[i]],
                          EOO = all.crit$EOO,
                          AOO = all.crit$AOO,
                          decline = all.crit$declineB,
                          pop.reduction = all.crit$reduction_A12,
                          pop.size = all.crit$pop.size,
                          pop.size.low = all.crit$pop.size.low, 
                          locations = all.crit$nbe_loc_total,
                          sever.frag = all.crit$sever.frag,
                          ext.fluct = NULL,
                          subpop = all.crit$Nbe_subPop, 
                          subcriteria = subcriteria[i])
}

rm(subcriteria)

##### CONSENSUS ASSESSMENT #####

assess.df <- all.crit[, c("species", "A1", "A2", "B1", "B2", "C1", "C2", "D")]
tmp <- ConR:::cat_mult_criteria(assess.df)
table(tmp$species == all.crit$species)
all.crit <- cbind.data.frame(all.crit, 
                             tmp[, c("category", "main.criteria", "aux.criteria")], 
                             stringsAsFactors = FALSE)
rm(assess.df)

##### Loading previous assessment using only inventory data #####
prev.assess <- read.csv('SI_assessment_sc.csv')
names(prev.assess)[37] <- "category_inv"

tmp <- merge(all.crit[,c("species", "category")],
             prev.assess[, c("species", "category_inv")],
             by = "species", all.x = TRUE, sort = FALSE)

tmp <- tmp[order(tmp$species),]
table(all.crit$species == tmp$species)

all.crit <- cbind.data.frame(all.crit, 
                             tmp[, c("category_inv")], 
                             stringsAsFactors = FALSE)

##### DOWNLISTING #####
rarity = read.csv('data_old_stuff/rarity_iffsc.csv', sep = ';')
rarity = rarity[, c(1:2)]
rarity = rarity[-c(130), ] #double Inga vera

#Correct names from rarity df
syn.br <- read.csv("new_synonyms_floraBR.csv", na.strings = c("", " ", NA), as.is = TRUE)
syn.br <- syn.br[syn.br$status %in% c("replace", "invert"), ]
for (i in 1:dim(syn.br)[1]) {
  sp.i <- syn.br$original.search[i]
  rpl.i <- syn.br$search.str[i]
  st.i <- syn.br$status[i]
  
  if (st.i == "replace")
    rarity$spp[rarity$spp %in% sp.i] <- rpl.i
}

#Correct names from rarity df
syn.br <- read.csv("new_synonyms_IFFSC.csv", sep = ';')
syn.br <- syn.br[syn.br$status %in% c("replace", "invert"), ]
for (i in 1:dim(syn.br)[1]) {
  sp.i <- syn.br$original.search[i]
  rpl.i <- syn.br$search.str[i]
  st.i <- syn.br$status[i]
  
  if (st.i == "replace")
    rarity$spp[rarity$spp %in% sp.i] <- rpl.i
}

#Since after correcting names we got double A. emarginata, we must remove one
rarity <- rarity[order(rarity$spp), ]
rarity <- rarity[-c(17), ]

#Other names are correct?
length(unique(rarity$spp)) #yes

all.crit2 = merge(all.crit, rarity, by.x = "species", by.y = "spp", all.x = T, sort = FALSE)
all.crit2$downlist[is.na(all.crit2$downlist)] <- "no"
table(all.crit2$downlist)
tmp = all.crit2$category[grepl("yes", all.crit2$downlist) & !all.crit2$category %in% c("NT", "LC")]
tmp1 = ConR:::cat_downlist(tmp)
all.crit2$category[grepl("yes", all.crit2$downlist) & !all.crit2$category %in% c("NT", "LC")] <- tmp1
table(all.crit2$category)

#Cleaning
all.crit2$category <- gsub("o$", "", all.crit2$category)
table(all.crit2$category)

#### SAVING THE FINAL RESULTS TABLE ####
all.crit2 <- all.crit2[order(all.crit2$species), ]
all.crit2 <- all.crit2[order(all.crit2$species), ]
#saveRDS(all.crit.new, "Resultados/all.criteria.cnbot.rds")
saveRDS(all.crit2, "data/all.criteria_herb.rds")
write.csv(all.crit2, "SI_assessment_sc_herb.csv")










