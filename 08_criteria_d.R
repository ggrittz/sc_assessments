rm(list=ls())

library(tidyverse)
library(ConR)
library(red)
library(circlize)
library(flora)

##### ASSESSING VERY SMALL POPULATIONS - CRITERIA D #####

#Already with population sizes estimated for all years
res.means <- readRDS("data/threat_mean_pop_sizes_infer.rds")
mean.pop.sizes <- readRDS("data/threat_mean_pop_sizes_for_ConR.rds")
low.pop.sizes <- readRDS("data/threat_low_pop_sizes_for_ConR.rds")
high.pop.sizes <- readRDS("data/threat_high_pop_sizes_for_ConR.rds")

#res.means with two species with wrong.names
#names(res.means)[[178]] <- "Erythrina crista-galli"
#names(res.means)[[583]] <- "Solanum sanctae-catharinae"
#saveRDS(res.means, "data/threat_mean_pop_sizes_infer.rds")

#Putting data in the ConR format
spp <- names(res.means)
nrows <- length(spp)
decline.models <- matrix(NA, ncol = 2, nrow = nrows,
                         dimnames = list(spp, c("Before_1985", "After_1985")))

for(x in 1:length(res.means)) {
  decline.models[x, 1] <- unique(res.means[[x]]$Modelo[res.means[[x]]$Year < 1985])
  decline.models[x, 2] <- unique(res.means[[x]]$Modelo[res.means[[x]]$Year > 1985])
}

##### LOADING HABITAT AND ECOLOGY DATA #####
#Includes species info on Generation Length and Proportion of mature individuals
hab <- read.csv("treeco/threat_habitats.csv", as.is = TRUE)

#Adjusting hab names
#Correcting new synonyms
syn.br <- read.csv("new_synonyms_floraBR.csv", na.strings = c("", " ", NA), as.is = TRUE)
syn.br <- syn.br[syn.br$status %in% c("replace", "invert"), ]
for (i in 1:dim(syn.br)[1]) {
  sp.i <- syn.br$original.search[i]
  rpl.i <- syn.br$search.str[i]
  st.i <- syn.br$status[i]
  
  if (st.i == "replace")
    hab$Name_submitted[hab$Name_submitted %in% sp.i] <- rpl.i
}

#Correcting new synonyms
syn.br <- read.csv("new_synonyms_IFFSC.csv", sep = ';')
syn.br <- syn.br[syn.br$status %in% c("replace", "invert"), ]
for (i in 1:dim(syn.br)[1]) {
  sp.i <- syn.br$original.search[i]
  rpl.i <- syn.br$search.str[i]
  st.i <- syn.br$status[i]
  
  if (st.i == "replace")
    hab$Name_submitted[hab$Name_submitted %in% sp.i] <- rpl.i
}

#Merging
PopData <- merge(decline.models, hab[, c(2, 6, 14)], 
                 by.x= "row.names", by.y = "Name_submitted", all.x = TRUE)
names(PopData)[1] <- "species"
PopData <- PopData[order(PopData$species), ]

#Adjusting some missing information different names from threat_habitats.csv or missing names
#If generation length and/or p.est are missing, calculate the median value from the genre
#Adding info for missing species on threat_habitats.csv

#Adding info for Erythrina crista gali (missing hifen on my data)
PopData[135, 1] <- "Erythrina crista-galli"
PopData[135, 4] <- 50
PopData[135, 5] <- 0.4498

#Adding missing info for Solanum sanctae catharinae (missing hifen)
PopData[506, 1] <- "Solanum sanctae-catharinae"
PopData[506, 4] <- 50
PopData[506, 5] <- 0.4498

#Missing info for Monteverdia schumanniana
#missing.info = hab %>% filter(stringr::str_detect(Name_submitted, "Vernonanthura"))

#Check NAs
sum(is.na(PopData))

#Correcting wrong names on my data
low.pop.sizes[453, 1] <- "Solanum sanctae-catharinae"
low.pop.sizes[129, 1] <- "Erythrina crista-galli"
mean.pop.sizes[453, 1] <- "Solanum sanctae-catharinae"
mean.pop.sizes[129, 1] <- "Erythrina crista-galli"
high.pop.sizes[453, 1] <- "Solanum sanctae-catharinae"
high.pop.sizes[129, 1] <- "Erythrina crista-galli"

#Verifying duplicate names
n_occur <- data.frame(table(PopData$species))
n_occur[n_occur$Freq > 1, ]
PopData[PopData$species %in% n_occur$Var1[n_occur$Freq > 1], ]

PopData <- aggregate(PopData[, 4:5], list(PopData$species), mean)
names(PopData)[1] <- "species"

#Filtering the populational datasets
ids <- mean.pop.sizes$species %in% PopData$species
mean.pop.sizes <- mean.pop.sizes[ids, ]
low.pop.sizes <- low.pop.sizes[ids, ]
table(mean.pop.sizes$species == PopData$species)
rm(res.means)

##### APPLYING CRITERIA D #####
#Getting previous assessments from Flora do Brasil and CONSEMA/2014

#Consema/2014
consema <- read.csv("consema_2014.csv", header = TRUE, sep = ';')

#Flora do Brasil
florinha <- get.taxa(unique(PopData$species))
PopData$species <- as.character(PopData$species)

previous_assess <- merge(PopData, consema, by = "species", all.x = TRUE, sort = FALSE)
previous_assess <- merge(previous_assess, florinha[, c("original.search", "threat.status")], by.x = "species",
                         by.y = "original.search", all.x = TRUE, sort = FALSE)

previous_assess <- previous_assess[, c(4, 1, 5, 6)]
names(previous_assess) <- c("family", "species", "consema", "flora")
saveRDS(previous_assess, 'data/previous_assessments_spp.rds')

##### Reading the files needed from species distributions #####
spp = readRDS('data/previous_assessments_spp.rds')

EOO <- readRDS("data/04_EOO_numeric.rds")
EOO <- tibble::rownames_to_column(EOO, "tax")
EOO$tax <- as.character(EOO$tax)
names(EOO)[1] <- "species"

AOO <- readRDS("data/04_AOO.rds")
AOO <- tibble::rownames_to_column(AOO, "tax")
AOO$tax <- as.character(AOO$tax)
names(AOO)[1] <- "species"

localities <- readRDS("data/04_number_of_locations.rds")

spp1<- merge(EOO, spp, by = "species", all.x = TRUE, sort = FALSE)
spp1<- merge(spp1, AOO, by = "species", all.x = TRUE, sort = FALSE)
spp1<- merge(spp1, localities, by = "species", all.x = TRUE, sort = FALSE)

#Creating and merging the df for analysis
df <- data.frame(species = mean.pop.sizes$species,
                 pop.size = mean.pop.sizes$`2019`,
                 pop.size.low = low.pop.sizes$`2019`,
                 p = PopData$p.est, stringsAsFactors = FALSE)

df1 <- merge(df, spp1, by = "species",
             all.x = TRUE, sort = FALSE)

df1 <- df1[order(df1$species), ]


##### ASSESSMENTS #####
#Optimal paramemeters
critD <- ConR::criterion_D(pop.size = df1$pop.size, 
                           Name_Sp = df1$species, 
                           AOO = df1$AOO,
                           n.Locs = df1$OutsideUC, #utilizando dados mais amplos para o criterio D
                           prop.mature = df1$p,
                           subcriteria = c("D", "D2"),
                           D.threshold = c(1000, 250, 50), 
                           AOO.threshold = 16, Loc.threshold = 2)

table(critD$D)
#saveRDS(critD, "data/critD.rds")

#Varying values of p
ps <- sort(c(1, .85, .72, .60, .49, .51, .58, .31, 
             .25, .33, .64, .45, .35, .28, 0.18, 0.4))
critD.all <- cbind.data.frame(critD,
                              D = as.character(criterion_D(df1$pop.size, 
                                                           Name_Sp = df1$species, 
                                                           subcriteria = c("D"),
                                                           prop.mature = ps[1])[, c("D")]), 
                              stringsAsFactors = FALSE)
for(i in 2:length(ps)){
  critD.all <- cbind.data.frame(critD.all,
                                D = as.character(criterion_D(df1$pop.size, 
                                                             Name_Sp = df1$species, 
                                                             subcriteria = c("D"),
                                                             prop.mature = ps[i])[, c("D")]), 
                                stringsAsFactors = FALSE)
}


##### Calculating the Red List Index for subcriterion A1 and A2 #####
#Optmimal parameters
all.GL2 <- critD.all[, c(1:4, 6:9, 5, 10:25)]
for(i in 9:25) all.GL2[, i] <- as.character(all.GL2[, i])
for(i in 9:25) all.GL2[, i] <- gsub("LC or NT", "LC", all.GL2[, i])

#RLI
rli.all2 <- apply(all.GL2[, 9:25], 2, red::rli, boot = TRUE, runs = 4999)
apply(all.GL2[, 9:25], 2, table)

#Renaming columns
names(all.GL2)[grepl("D\\.[0-9]", names(all.GL2))] <- paste0("D.p", ps, sep = "")

## Renaming the LC category
all.GL2[] <- lapply(all.GL2, gsub, pattern = "^LC$", replacement = "LC or NT")

## Adding the low population estimates
table(all.GL2$species == df1$species)
all.GL2$pop.size.low <- df1$pop.size.low * df1$p

#Saving
saveRDS(all.GL2, "data/criterionD_all_prop_mature.rds")


#### FIGURES ####
jpeg(filename = "figures/Figure_SZ.jpg", width = 2500, height = 2000, units = "px", pointsize = 12,
     res = 300, family = "sans", type = "cairo", bg = "white", quality = 100)
#par(mfrow = c(1, 2))
par(mar=c(4, 4, 0.75, 0.5), mgp = c(2.5, 0.25, 0),tcl = -0.2,las = 1)
#Optimum generation length
plot(rev(rli.all2[2, grepl("D\\.[0-9]", colnames(rli.all2))]) ~ rev(ps), #type = "b",
     #xaxp = c(1890, 2018, 5), yaxp = c(60, 220, 6), 
     xaxt = "n",# yaxt = "n", 
     cex.lab = 1.2,
     xlab = "Prop. mature individuals", ylab = "Red List Index", 
     pch = 19, ylim = c(0.99, 1))
#axis(1, at = rev(ps), cex.axis = 1)
axis(1, at=seq(0.2, 1, 0.1), cex.axis = 1)
#axis(2, at=c(60, 80, 100, 120, 140, 160, 180, 200, 220), cex.axis = 1)
arrows(x0=rev(ps), y0 = rev(rli.all2[1,grepl("D\\.[0-9]", colnames(rli.all2))]), #using mean CIs
       y1 = rev(rli.all2[3,grepl("D\\.[0-9]", colnames(rli.all2))]),
       code = 3, angle = 90, length = 0.05)
#arrows(x0=rev(ps), y0 = rev(rli.all1[1,grepl("\\.p", colnames(rli.all1))][17:32]), #using low and high CIs
#       y1 = rev(rli.all1[3,grepl("\\.p", colnames(rli.all1))][33:48]),
#       code = 3, angle = 90, length = 0.05, col = 2)
#legend("topleft", expression(bold(A)), bty="n", cex = 1.3)
abline(h = rli.all2[2,1], lty = 2)
legend("bottomright", c("Group-specific", "Fixed"),
       lty = c(3, 0), pch = c(NA, 19),
       bty = "n", lwd = 2)
dev.off()

