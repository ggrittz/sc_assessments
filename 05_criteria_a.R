rm(list=ls())

##### CRITERIA A: POPULATION DECLINE #####

library(ConR)
library(red)
library(circlize)
library(flora)

##### LOADING THREAT POPULATION SIZE DATA (TREECO) #####
#Already with pop. sizes estimated for all necessary years
res.means <- readRDS("data/threat_mean_pop_sizes_infer.rds")
low.pop.sizes <- readRDS("data/threat_low_pop_sizes_for_ConR.rds")
mean.pop.sizes <- readRDS("data/threat_mean_pop_sizes_for_ConR.rds")
high.pop.sizes <- readRDS("data/threat_high_pop_sizes_for_ConR.rds")

#Putting data in the ConR format
spp <- names(res.means)
nrows <- length(spp)
decline.models <- matrix(NA, ncol = 2, nrow = nrows,
                         dimnames = list(spp, c("Before 1985", "After 1985")))

for(x in 1:length(res.means)) {
  decline.models[x, 1] <- unique(res.means[[x]]$Modelo[res.means[[x]]$Year < 1985])
  decline.models[x, 2] <- unique(res.means[[x]]$Modelo[res.means[[x]]$Year > 1985])
}

##### LOADING THREAT HABITAT AND ECOLOGY DATA #####
#Includes species info on Generation Length and Proportion of mature individuals
hab <- read.csv("treeco/threat_habitats.csv", as.is = TRUE)

#Adjusting hab names
#Correcting new synonyms
syn.br <- read.csv("new_synonyms_floraBR.csv", na.strings = c(""," ",NA), as.is = TRUE)
syn.br <- syn.br[syn.br$status %in% c("replace", "invert"), ]
for (i in 1:dim(syn.br)[1]) {
  sp.i <- syn.br$original.search[i]
  rpl.i <- syn.br$search.str[i]
  st.i <- syn.br$status[i]
  
  if (st.i == "replace")
    hab$Name_submitted[hab$Name_submitted %in% sp.i] <- rpl.i
}

#Now IFFSC names
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
PopData <- merge(decline.models, hab[, c(2, 6, 14)], by.x= "row.names", by.y = "Name_submitted", all.x = TRUE)
names(PopData)[1] <- "species"
PopData <- PopData[order(PopData$species), ]

#Synonyms created more than 1 new value with the same name?
n_occur <- data.frame(table(PopData$species))
dup.freq = n_occur[n_occur$Freq > 1, ]

#Adding info for Erythrina crista gali (missing hifen on my data)
PopData[135, 1] <- "Erythrina crista-galli"
PopData[135, 4] <- 50
PopData[135, 5] <- 0.4498

#Adding missing info for Solanum sanctae catharinae (missing hifen)
PopData[506, 1] <- "Solanum sanctae-catharinae"
PopData[506, 4] <- 50
PopData[506, 5] <- 0.4498

#Take mean from generation length and p.est when multiple names occurred because of synonymization
PopData <- aggregate(PopData[, 4:5], list(PopData$species), mean)
names(PopData)[1] <- "species"

PopData <- PopData[order(PopData$species), ]

## Filtering the populational datasets
#Correcting those names again...
low.pop.sizes[129, 1] <- "Erythrina crista-galli"
mean.pop.sizes[129, 1] <- "Erythrina crista-galli"
high.pop.sizes[129, 1] <- "Erythrina crista-galli"

low.pop.sizes[453, 1] <- "Solanum sanctae-catharinae"
mean.pop.sizes[453, 1] <- "Solanum sanctae-catharinae"
high.pop.sizes[453, 1] <- "Solanum sanctae-catharinae"

ids <- mean.pop.sizes$species %in% PopData$species
mean.pop.sizes <- mean.pop.sizes[ids, ]
low.pop.sizes <- low.pop.sizes[ids, ]
high.pop.sizes <- high.pop.sizes[ids, ]

table(mean.pop.sizes$species == PopData$species)
table(low.pop.sizes$species == PopData$species)
table(high.pop.sizes$species == PopData$species)

rm(res.means)

PopData$GenerationLength.range <- round(PopData$GenerationLength.range, digits = 0)
PopData$GenerationLength.range <- as.integer(PopData$GenerationLength.range)

#### APPLYING CRITERIA A ####
## ASSESSMENTS USING SPECIES-SPECIFIC PROXIES OF GENERATION LENGTH ##
critA <- criterion_A(mean.pop.sizes, assess.year = 2019,
                     project.years = NULL, subcriteria = c("A1", "A2"),
                     generation.time = PopData$GenerationLength.range)

critA.low <- criterion_A(low.pop.sizes, assess.year = 2019,
                         project.years = NULL, subcriteria = c("A1","A2"),
                         generation.time = PopData$GenerationLength.range)

critA.high <- criterion_A(high.pop.sizes, assess.year = 2019,
                         project.years = NULL, subcriteria = c("A1","A2"),
                         generation.time = PopData$GenerationLength.range)

table(critA$A1)
table(critA.low$A1)
table(critA.high$A1)

saveRDS(critA.low, 'data/critA_low.rds')
## Calculating the Red List Index for subcriterion A1 and A2
#rli.all <- apply(critA, 2, red::rli, boot = TRUE, runs = 4999)


## ASSESSMENTS of the influence of different GENERATION LENGTHS (FIXED FOR ALL SPECIES) ##
all.GL <- cbind.data.frame(critA,
                           criterion_A(mean.pop.sizes, assess.year = 2019, subcriteria = c("A1","A2"),
                                       generation.time = 10)[,c("reduction_A12","A1","A2")],
                           criterion_A(mean.pop.sizes, assess.year = 2019, subcriteria = c("A1","A2"),
                                       generation.time = 20)[,c("reduction_A12","A1","A2")],
                           criterion_A(mean.pop.sizes, assess.year = 2019, subcriteria = c("A1","A2"),
                                       generation.time = 25)[,c("reduction_A12","A1","A2")],
                           criterion_A(mean.pop.sizes, assess.year = 2019, subcriteria = c("A1","A2"),
                                       generation.time = 30)[,c("reduction_A12","A1","A2")],
                           criterion_A(mean.pop.sizes, assess.year = 2019, subcriteria = c("A1","A2"),
                                       generation.time = 35)[,c("reduction_A12","A1","A2")],
                           criterion_A(mean.pop.sizes, assess.year = 2019, subcriteria = c("A1","A2"),
                                       generation.time = 40)[,c("reduction_A12","A1","A2")],
                           criterion_A(mean.pop.sizes, assess.year = 2019, subcriteria = c("A1","A2"),
                                       generation.time = 45)[,c("reduction_A12","A1","A2")],
                           criterion_A(mean.pop.sizes, assess.year = 2019, subcriteria = c("A1","A2"),
                                       generation.time = 50)[,c("reduction_A12","A1","A2")],
                           criterion_A(mean.pop.sizes, assess.year = 2019, subcriteria = c("A1","A2"),
                                       generation.time = 55)[,c("reduction_A12","A1","A2")],
                           stringsAsFactors = FALSE, deparse.level = 0
)

all.GL <- all.GL[,c(1:4,6,7,10,
                    5,11,14,17,20,23,26,29,32,35, #Reductions A12
                    8,12,15,18,21,24,27,30,33,36, #A1
                    9,13,16,19,22,25,28,31,34,37)] #A2
for(i in 18:37) all.GL[,i] <- as.character(all.GL[,i])
for(i in 18:37) all.GL[,i] <- gsub("LC or NT", "LC", all.GL[,i])

## Calculating the Red List Index for subcriterion A1 and A2
rli.all <- apply(all.GL[,19:37], 2, red::rli, boot = TRUE, runs = 4999)

###################
#### FIGURE SX ####
###################

jpeg(filename = "figures/Figure_SX.jpg", width = 4000, height = 2000, units = "px", pointsize = 12,
     res = 300, family = "sans", type = "cairo", bg = "white", quality = 100)
gls = c(10,20,25,30,35,40,45,50,55)
par(mfrow=c(1,2))
par(mar=c(3,3.5,0.75,0.5), mgp=c(1.9,0.25,0),tcl=-0.2,las=1)
#subcriterion A1
plot(rli.all[2,grepl("A1\\.", colnames(rli.all))] ~ gls, #type = "b",
     #xaxp = c(1890, 2018, 5), yaxp = c(60, 220, 6), 
     xaxt = "n",# yaxt = "n", 
     cex.lab = 1.2,
     xlab = "Generation lenght (years)", ylab = "Red List Index", 
     pch=19, ylim = c(0.39,1))
axis(1, at = gls, cex.axis = 1)
#axis(2, at = c(60,80,100,120,140,160,180,200,220), cex.axis = 1)
arrows(x0 = gls, y0 = rli.all[1,grepl("A1\\.", colnames(rli.all))], 
       y1 = rli.all[3,grepl("A1\\.", colnames(rli.all))],
       code = 3, angle = 90, length = 0.05)
legend("topright", expression(bold(A)), bty="n", cex = 1.3)
abline(h=rli.all[2,1], lty = 2)

#Subcriterion A2
plot(rli.all[2,grepl("A2\\.", colnames(rli.all))] ~ gls, #type = "b",
     #xaxp = c(1890, 2018, 5), yaxp = c(60, 220, 6), 
     xaxt = "n", #yaxt = "n", 
     cex.lab = 1.2,
     xlab = "Generation lenght (years)", ylab = "Red List Index", 
     pch=19, ylim = c(0.39,1))
axis(1, at = gls, cex.axis = 1)
arrows(x0 = gls, y0 = rli.all[1,grepl("A2\\.", colnames(rli.all))], 
       y1 = rli.all[3,grepl("A2\\.", colnames(rli.all))],
       code = 3, angle = 90, length = 0.05)
legend("topright", expression(bold(B)), bty = "n", cex = 1.3)
abline(h = rli.all[2,11], lty = 2)
legend("bottomleft", c("Group-specific", "Fixed"),
       lty = c(2,0), pch = c(NA,19),
       bty = "n", lwd = 2)
dev.off()

#Where are the changes related to low generation lengths (25 years)
critA.gl25 <- criterion_A(mean.pop.sizes,
                          assess.year = 2019,
                          project.years = NULL,
                          subcriteria = c("A1","A2"),
                          generation.time = 25)

critA.gl25[critA$A1 %in% "CR" & critA.gl25$A1 %in% "EN", ] #0 cases
critA.gl25[critA$A1 %in% "CR" & critA.gl25$A1 %in% "VU", ] #0 cases
critA.gl25[critA$A1 %in% "CR" & critA.gl25$A1 %in% "LC or NT", ] #0 cases
critA.gl25[critA$A1 %in% "EN" & critA.gl25$A1 %in% "VU", ] #1 case
critA.gl25[critA$A1 %in% "EN" & critA.gl25$A1 %in% "LC or NT", ] #11 cases
critA.gl25[critA$A1 %in% "VU" & critA.gl25$A1 %in% "LC or NT", ] #331 cases

#Inpecting some classic examples
critA[critA$species %in% "Araucaria angustifolia", ]
critA.gl25[critA$species %in% "Araucaria angustifolia", ]
critA[critA$species %in% "Euterpe edulis", ]
critA.gl25[critA$species %in% "Euterpe edulis", ]
critA[critA$species %in% "Cedrela fissilis",]
critA.gl25[critA$species %in% "Cedrela fissilis", ]
critA[critA$species %in% "Cecropia pachystachya", ]
critA.gl25[critA$species %in% "Cecropia pachystachya", ]

#Renaming columns
names(all.GL)[grepl("\\.1$", names(all.GL))] <- gsub("\\.1$", ".10ys", names(all.GL)[grepl("\\.1$", names(all.GL))])
names(all.GL)[grepl("\\.2$", names(all.GL))] <- gsub("\\.2$", ".20ys", names(all.GL)[grepl("\\.2$", names(all.GL))])
names(all.GL)[grepl("\\.3$", names(all.GL))] <- gsub("\\.3$", ".25ys", names(all.GL)[grepl("\\.3$", names(all.GL))])
names(all.GL)[grepl("\\.4$", names(all.GL))] <- gsub("\\.4$", ".30ys", names(all.GL)[grepl("\\.4$", names(all.GL))])
names(all.GL)[grepl("\\.5$", names(all.GL))] <- gsub("\\.5$", ".35ys", names(all.GL)[grepl("\\.5$", names(all.GL))])
names(all.GL)[grepl("\\.6$", names(all.GL))] <- gsub("\\.6$", ".40ys", names(all.GL)[grepl("\\.6$", names(all.GL))])
names(all.GL)[grepl("\\.7$", names(all.GL))] <- gsub("\\.7$", ".45ys", names(all.GL)[grepl("\\.7$", names(all.GL))])
names(all.GL)[grepl("\\.8$", names(all.GL))] <- gsub("\\.8$", ".50ys", names(all.GL)[grepl("\\.8$", names(all.GL))])
names(all.GL)[grepl("\\.9$", names(all.GL))] <- gsub("\\.9$", ".55ys", names(all.GL)[grepl("\\.9$", names(all.GL))])

#Renaming the LC category
all.GL[] <- lapply(all.GL, gsub, pattern = "^LC$", replacement = "LC or NT")

##### Saving #####
saveRDS(all.GL, "data/criterionA_all_GLs.rds")

rm(list=ls())
