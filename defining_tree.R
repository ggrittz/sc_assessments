library(tidyverse)
library(plantR)
library(flora)

#First filters
florasc_gbif = read_tsv('C:/Users/Master/OneDrive - FURB/Mestrado/data/florasc_gbif.csv')
florasc_gbif = florasc_gbif %>% select(class, family, species)
florasc_gbif = na.omit(florasc_gbif)

#Selecting only unique rows
florasc = distinct(florasc_gbif, species, .keep_all = TRUE)

#Correcting names, if needed
spp = unique(florasc$species)
#Flora do Brasil
fb = get.taxa(spp, life.form = TRUE)
fb = fb %>% select(family, search.str, notes, life.form)
fb = na.omit(fb)
fb = fb %>% select(family, search.str, life.form)
names(fb) = c("family", "species", "life_form")

trees_arbs = filter(fb, grepl('Árvore', life_form))
table(trees_arbs$life_form)

#Removing duplicates
trees_arbs = trees_arbs %>%
  group_by(species) %>%
  slice(1)

###Loading IFFSC data to compare###
#Regeneration plots already have the corrected names
iffsc_reg = read.csv('C:/Users/Master/OneDrive - FURB/Mestrado/iffsc_regen_ago21.csv', header = TRUE, sep = ';')
iffsc_reg$spp = word(iffsc_reg$Espécie, start = 1, end = 2)
iffsc_reg = na.omit(iffsc_reg)
#Number of times the species occur in regeneration
n_reg = table(iffsc_reg$spp)
n_reg = as.data.frame(n_reg)
n_reg$Var1 = as.character(n_reg$Var1)


#Arboreal data needs name correction
iffsc_arb = read.csv('C:/Users/Master/OneDrive - FURB/Mestrado/iffsc_geral_jan18.csv', header = TRUE, sep = ';')
iffsc_arb$Espécie = word(iffsc_arb$Espécie, start = 1, end = 2)
iffsc_arb = na.omit(iffsc_arb)
#Number of times the species occur in arboreal
n_arb = table(iffsc_arb$Espécie)
n_arb = as.data.frame(n_arb)
n_arb$Var1 = as.character(n_arb$Var1)
#Correcting names
fb = get.taxa(n_arb$Var1)
n_arb$Var1 = fb$search.str

trees_list = merge(trees_arbs, n_reg, by.x = "species", by.y = "Var1", all.x = TRUE, sort = F)
names(trees_list)[4] <- "n_reg"
trees_list = merge(trees_list, n_arb, by.x = "species", by.y = "Var1", all.x = TRUE, sort = F)
names(trees_list)[5] <- "n_arb"

#write.csv(trees_list, 'C:/Users/Master/OneDrive - FURB/Mestrado/sc_trees_list.csv')

#Removing exotics
treelist = read.csv('C:/Users/Master/OneDrive - FURB/Mestrado/sc_trees_list.csv', encoding = 'UTF-8-BOM', sep = ';')
treelist = treelist[, -c(1)]
length(unique(treelist$species))

#Removing duplicate terms (if any)
n_occur <- data.frame(table(treelist$species))
n_occur[n_occur$Freq > 1,]
treelist[treelist$species %in% n_occur$Var1[n_occur$Freq > 1],]

#Getting establishment...
floralist = get.taxa(treelist$species, establishment = TRUE)

final = merge(treelist, floralist[, c("original.search", "establishment")], 
              by.x = "species", 
              by.y = "original.search",
              all.x = TRUE, sort = F)
final = filter(final, grepl("native", establishment))

write_csv(final, 'C:/Users/Master/OneDrive - FURB/Mestrado/sc_trees_list.csv')

#### NEW STUFF ####
sc_treelist = read.csv('C:/Users/Master/OneDrive - FURB/Mestrado/sc_trees_list.csv', fileEncoding = 'UTF-8-BOM')
sc_reviewed = read.csv2('C:/Users/Master/OneDrive - FURB/Mestrado/Para revisão/reviewed_cr_spp.csv')

sc_treelist_notes = merge(sc_treelist, sc_reviewed[, c("mainly.as.a.shrub.treelet", "dado", "ref", "species")],
                          by.x = "species", by.y = "species", all.x = TRUE, sort = F)


#Does it occur in SC in Flora2020?
spp_flora = get.taxa(spp, states = T)
spp_flora_sc = filter(spp_flora, grepl("SC", occurrence))

sc_treelist_notes = merge(sc_treelist_notes, spp_flora_sc[, c("original.search", "occurrence")], by.x = "species",
                          by.y = "original.search", all.x = TRUE, sort = F)

#Ocorre ou não?
sc_treelist_notes$occurrence = sc_treelist_notes$occurrence %>% replace_na("não ocorre")
sc_treelist_notes$occurrence[ sc_treelist_notes$occurrence != "não ocorre" ] <- "ocorre"

#Quantity of occurrences in GBIF
gbif = read_csv2('C:/Users/Master/Desktop/0013562-210819072339941/gbif_again.csv')
gbif$abundance = 1
gbif = filter(gbif, taxonRank == "SPECIES")
abs = gbif %>% group_by(species, abundance) %>%
  summarize(ab = sum(abundance))
abs = abs[, c("species", "ab")]

#merging again...
sc_treelist_notes = merge(sc_treelist_notes, abs, by = "species", all.x = TRUE, sort = F)
sc_treelist_final = sc_treelist_notes[, -c(11)]

write.csv(sc_treelist_final, 'C:/Users/Master/OneDrive - FURB/Mestrado/sc_trees_list_again.csv')
