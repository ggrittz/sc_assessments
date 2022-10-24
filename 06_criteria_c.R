##### EVALUATING POPULATION DECLINE (CRITERION C) #####
rm(list=ls())

library(dplyr)
library(ConR)
library(red)
library(circlize)
library(flora)
library(plyr)

##### LOADING POPULATION SIZE DATA (SC) #####
#Already with population sizes estimated for all years
res.means <- readRDS("data/threat_mean_pop_sizes_infer.rds")
mean.pop.sizes <- readRDS("data/threat_mean_pop_sizes_for_ConR.rds")
low.pop.sizes <- readRDS("data/threat_low_pop_sizes_for_ConR.rds")
high.pop.sizes <- readRDS("data/threat_high_pop_sizes_for_conR.rds")

#Putting data in the ConR format
spp <- names(res.means)
nrows <- length(spp)
decline.models <- matrix(NA, ncol = 2, nrow = nrows,
                         dimnames = list(spp, c("Before_1985", "After_1985")))
for(x in 1:length(res.means)) {
  decline.models[x, 1] <- unique(res.means[[x]]$Modelo[res.means[[x]]$Year<1985])
  decline.models[x, 2] <- unique(res.means[[x]]$Modelo[res.means[[x]]$Year>1985])
}

##### LOADING HABITAT AND ECOLOGY DATA #####
#Includes species info on generation length and proportion of mature individuals
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

#Correcting new synonyms (Now IFFSC)
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

#Mean from multiple equal names
PopData <- aggregate(PopData[, 4:5], list(PopData$species), mean)
names(PopData)[1] <- "species"


#Filtering the populational datasets
ids <- mean.pop.sizes$species %in% PopData$species
mean.pop.sizes <- mean.pop.sizes[ids, ]
low.pop.sizes <- low.pop.sizes[ids, ]
high.pop.sizes <- high.pop.sizes[ids, ]

table(low.pop.sizes$species == PopData$species)
table(mean.pop.sizes$species == PopData$species)
table(high.pop.sizes$species == PopData$species)

rm(res.means)

#Saving corrected names
saveRDS(low.pop.sizes, "data/threat_low_pop_sizes_for_ConR.rds")
saveRDS(mean.pop.sizes, "data/threat_mean_pop_sizes_for_ConR.rds")
saveRDS(high.pop.sizes, "data/threat_high_pop_sizes_for_ConR.rds")

##### CRITERIA C #####
# Getting info on the number of subpopulations for each species
subpops <- readRDS("data/04_n_subpop.rds")
subpops <- tibble::rownames_to_column(subpops, "species")
subpops <- subpops[order(subpops$species), ]

#Final check
table(subpops$species == PopData$species)
table(low.pop.sizes$species == subpops$species)

#Doing the same for mean.pop.sizes and low.pop.sizes
subpops <- merge(low.pop.sizes[,c("species","2019")], subpops,
                 by = "species", all.x = TRUE, sort = FALSE)
subpops <- subpops[order(subpops$species), ]
table(low.pop.sizes$species == subpops$species)
subpop.sizes <- vector("list", dim(subpops)[1])
names(subpop.sizes) <- low.pop.sizes$species

#Pop size in each subpopulation
for(i in 1:length(subpop.sizes)) {
  subpop.sizes[[i]] <- rep(subpops$`2019`[i]/subpops$Number_subpop[i],
                           subpops$Number_subpop[i])
}

#Testing
table(round(low.pop.sizes[, which(names(low.pop.sizes) == 2019)], 0) == round(sapply(subpop.sizes, sum), 0))

##### Running the criterion C for the optimal parameters (GL and p) #####
#Subsetting two species
#test.df.popsize = subset(low.pop.sizes, species %in% c("Araucaria angustifolia", "Aspidosperma australe"))
#test.df.hab = subset(hab, Name_submitted %in% c("Araucaria angustifolia", "Aspidosperma australe"))
#test.subpop.size = subpop.sizes[c("Araucaria angustifolia", "Aspidosperma australe")]

#ignore years not working! avisar Renato
critC_low <- criterion_C(x = low.pop.sizes[, c(1, 14:72)],
                         assess.year = 2019,
                         project.years = NULL,
                         project = FALSE,
                         #ignore.years = c(1719, 1734, 1749, 1764, 1779, 1794, 1809, 1819, 1824, 1829, 1839, 1849), #1850 upper limit
                         recent.year = 2000,
                         subcriteria = c("C1","C2"),
                         models = c("linear", "quadratic", "exponential", "logistic", "general_logistic"),
                         generation.time = PopData$GenerationLength.range,
                         prop.mature = PopData$p.est,
                         subpop.size = subpop.sizes,
                         parallel = TRUE,
                         NbeCores = 4)

#doing the same for mean.pop.sizes and low.pop.sizes
#subpops <- merge(mean.pop.sizes[,c("species","2019")], subpops,
#                 by = "species", all.x = TRUE, sort = FALSE)
#subpops <- subpops[order(subpops$species), ]
#table(mean.pop.sizes$species == subpops$species)
#subpop.sizes <- vector("list", dim(subpops)[1])
#names(subpop.sizes) <- mean.pop.sizes$species

#Pop size in each subpopulation
#for(i in 1:length(subpop.sizes)) {
# subpop.sizes[[i]] <- rep(subpops$`2019`[i]/subpops$Number_subpop[i],
#                           subpops$Number_subpop[i])
#}

#Testing
#table(round(mean.pop.sizes[, which(names(mean.pop.sizes) == 2019)], 0) == round(sapply(subpop.sizes, sum), 0))

#critC_mean <- criterion_C(x = mean.pop.sizes[, c(1, 14:72)],
#                     assess.year = 2019,
#                     project.years = NULL,
#                     project = FALSE,
#                     #ignore.years = c(1719, 1734, 1749, 1764, 1779, 1794, 1809, 1819, 1824, 1829, 1839, 1849), #1850 upper limit
#                     recent.year = 2000,
#                     subcriteria = c("C1", "C2"),
#                     models = c("linear", "quadratic", "exponential", "logistic", "general_logistic"),
#                     generation.time = PopData$GenerationLength.range,
#                     prop.mature = PopData$p.est,
#                     subpop.size = subpop.sizes,
#                     parallel = TRUE,
#                     NbeCores = 4)

#saveRDS(critC_mean, "data/critC_mean.rds")
saveRDS(critC_low, "data/critC_low.rds")


##### Running the criterion C for GL fixed at 25 years and the optimal p.est #####
critC.gl25 <- criterion_C(x = low.pop.sizes[, c(1, 14:72)],
                          assess.year = 2019,
                          project.years = NULL,
                          project = FALSE,
                          #ignore.years = c(1719, 1734, 1749, 1764, 1779, 1794, 1809, 1819, 1824, 1829, 1839, 1849),
                          recent.year = 2000,
                          subcriteria = c("C1","C2"),
                          models = c("linear", "quadratic", "exponential", "logistic", "general_logistic"),
                          generation.time = 25,
                          prop.mature = PopData$p.est,
                          subpop.size = subpop.sizes,
                          parallel = TRUE,
                          NbeCores = 4)

saveRDS(critC.gl25, "data/critC_GL_25ys_low.rds")
critC.gl25 <- readRDS("data/critC_GL_25ys_low.rds")

##### RUNNING ASSESSMENTS FOR DIFFERENT VALUES OF p #####
ps <- sort(c(1, .85, .72, .60, .49, .51, .58, .31, .25,
             .33, .64, .45, .35, .28, 0.18, 0.4), decreasing = TRUE)

#Optimal parameters
df <- cbind.data.frame(critC_low[,c("assess.pop.size", "cont.decline")],
                       critC_low[,grepl("reduction", names(critC_low))],
                       stringsAsFactors = FALSE)
df$assess.pop.size <- mean.pop.sizes$`2019`
df$assess.pop.size.low <- low.pop.sizes$`2019`
df$assess.pop.size.high <- high.pop.sizes$`2019`

res <- res.low <- res.high <- vector("list",length(ps))
names(res) <- names(res.low) <- names(res.high) <- ps
for(i in 1:length(ps)){
  C1 <- df
  C1$assess.pop.size <- df$assess.pop.size * ps[i]
  all_ranks <- cat_criterion_c(C1_df = C1, C2_df = NULL)
  res[[i]] <- as.character(all_ranks$ranks_C)
  
  C1$assess.pop.size <- df$assess.pop.size.low * ps[i]
  all_ranks <- cat_criterion_c(C1_df = C1, C2_df = NULL)
  res.low[[i]] <- as.character(all_ranks$ranks_C)
  
  C1$assess.pop.size <- df$assess.pop.size.high * ps[i]
  all_ranks <- cat_criterion_c(C1_df = C1, C2_df = NULL)
  res.high[[i]] <- as.character(all_ranks$ranks_C)
}
res <- do.call(cbind.data.frame, res)
colnames(res) <- paste0("C1.p",colnames(res))
res.low <- do.call(cbind.data.frame, res.low)
colnames(res.low) <- paste0("C1.p",colnames(res.low), ".low")
res.high <- do.call(cbind.data.frame, res.high)
colnames(res.high) <- paste0("C1.p",colnames(res.high), ".high")
critC.all <- cbind.data.frame(critC_low, res, res.low, res.high)

#Generation length = 25 years
df <- cbind.data.frame(critC.gl25[,c("assess.pop.size", "cont.decline")],
                       critC.gl25[,grepl("reduction", names(critC.gl25))],
                       stringsAsFactors = FALSE)
df$assess.pop.size <- mean.pop.sizes$`2019`
df$assess.pop.size.low <- low.pop.sizes$`2019`
df$assess.pop.size.high <- high.pop.sizes$`2019`

res <- res.low <- res.high <- vector("list",length(ps))
names(res) <- names(res.low) <- names(res.high) <- ps
for(i in 1:length(ps)){
  C1 <- df
  C1$assess.pop.size <- df$assess.pop.size * ps[i]
  all_ranks <- cat_criterion_c(C1_df = C1, C2_df = NULL)
  res[[i]] <- as.character(all_ranks$ranks_C)
  
  C1$assess.pop.size <- df$assess.pop.size.low * ps[i]
  all_ranks <- cat_criterion_c(C1_df = C1, C2_df = NULL)
  res.low[[i]] <- as.character(all_ranks$ranks_C)
  
  C1$assess.pop.size <- df$assess.pop.size.high * ps[i]
  all_ranks <- cat_criterion_c(C1_df = C1, C2_df = NULL)
  res.high[[i]] <- as.character(all_ranks$ranks_C)
}

res <- do.call(cbind.data.frame, res)
colnames(res) <- paste0("C1.p", colnames(res))
res.low <- do.call(cbind.data.frame, res.low)
colnames(res.low) <- paste0("C1.p", colnames(res.low), ".low")
res.high <- do.call(cbind.data.frame, res.high)
colnames(res.high) <- paste0("C1.p", colnames(res.high), ".high")
critC.all.gl25 <- cbind.data.frame(critC.gl25, res, res.low, res.high)

##### Calculating the Red List Index for subcriterion A1 and A2 #####
#Optimal parameters
library(red)

all.GL1 <- critC.all[, c(1:14, 15:62)]
for(i in 20:62) all.GL1[,i] <- as.character(all.GL1[,i])
for(i in 20:62) all.GL1[,i] <- gsub("LC or NT", "LC", all.GL1[,i])

rli.all1 <- apply(all.GL1[, 20:62], 2, red::rli, boot = TRUE, runs = 4999)

apply(all.GL1[, 20:62], 2, table)[1:2]

#Generation length = 25 years
all.GL1.gl25 <- critC.all.gl25[, c(1:14, 15:62)]
for(i in 20:62) all.GL1.gl25[, i] <- as.character(all.GL1.gl25[,i])
for(i in 20:62) all.GL1.gl25[, i] <- gsub("LC or NT", "LC", all.GL1.gl25[, i])

rli.all1.gl25 <- apply(all.GL1.gl25[, 20:62], 2, red::rli, boot = TRUE, runs = 4999)
apply(all.GL1.gl25[, 20:62], 2, table)[1:2]

#Renaming the LC category
all.GL1[] <- lapply(all.GL1, gsub, pattern = "^LC$", replacement = "LC or NT")
all.GL1.gl25[] <- lapply(all.GL1.gl25, gsub, pattern = "^LC$", replacement = "LC or NT")

#### Saving ####
saveRDS(all.GL1, "data/criterionC_all_prop_mature.rds")
saveRDS(all.GL1.gl25, "data/criterionC_all_prop_mature_gl25.rds")


###################
#### FIGURE SY ####
###################

jpeg(filename = "figures/Figure_SY.jpg", width = 4000, height = 2000, units = "px", pointsize = 12,
     res = 300, family = "sans", type = "cairo", bg = "white", quality = 100)
par(mfrow = c(1, 2))
par(mar = c(3, 3.5, 0.75, 0.5), mgp = c(1.9, 0.25, 0),tcl = -0.2, las = 1)
#Optimum generation length
plot(rev(rli.all1[2, grepl("\\.p", colnames(rli.all1))][1:16]) ~ rev(ps), #type = "b",
     #xaxp = c(1890, 2018, 5), yaxp = c(60, 220, 6), 
     xaxt = "n",# yaxt = "n", 
     cex.lab = 1.2,
     xlab = "Prop. mature individuals", ylab = "Red List Index", 
     pch = 19, ylim = c(0.95, 1))
#axis(1, at = rev(ps), cex.axis = 1)
axis(1, at = seq(0.2, 1, 0.1), cex.axis = 1)
#axis(2, at = c(60, 80, 100, 120, 140, 160, 180, 200, 220), cex.axis = 1)
arrows(x0 = rev(ps), y0 = rev(rli.all1[1,grepl("\\.p", colnames(rli.all1))][1:16]), #using mean CIs
       y1 = rev(rli.all1[3,grepl("\\.p", colnames(rli.all1))][1:16]),
       code = 3, angle = 90, length = 0.05)
#arrows(x0 = rev(ps), y0 = rev(rli.all1[1,grepl("\\.p", colnames(rli.all1))][17:32]), #using low and high CIs
#       y1 = rev(rli.all1[3,grepl("\\.p", colnames(rli.all1))][33:48]),
#       code = 3, angle = 90, length = 0.05, col = 2)
legend("bottomleft", expression(bold(A)), bty = "n", cex = 1.3)
abline(h = rli.all1[2,1], lty = 2)

#Generation length = 25 years
plot(rev(rli.all1.gl25[2, grepl("\\.p", colnames(rli.all1.gl25))][1:16]) ~ rev(ps), #type = "b",
     #xaxp = c(1890, 2018, 5), yaxp = c(60, 220, 6), 
     xaxt = "n",# yaxt = "n", 
     cex.lab = 1.2,
     xlab = "Prop. mature individuals", ylab = "Red List Index", 
     pch=19, ylim = c(0.95, 1))
#axis(1, at = rev(ps), cex.axis = 1)
axis(1, at = seq(0.2,1,0.1), cex.axis = 1)
#axis(2, at = c(60, 80, 100, 120, 140, 160, 180, 200, 220), cex.axis = 1)
arrows(x0 = rev(ps), y0 = rev(rli.all1.gl25[1,grepl("\\.p", colnames(rli.all1.gl25))][1:16]), #using mean CIs
       y1 = rev(rli.all1[3, grepl("\\.p", colnames(rli.all1.gl25))][1:16]),
       code = 3, angle = 90, length = 0.05)
# arrows(x0 = rev(ps), y0 = rev(rli.all1.gl25[1,grepl("\\.p", colnames(rli.all1.gl25))][17:32]), #using low and high CIs
#        y1 = rev(rli.all1.gl25[3, grepl("\\.p", colnames(rli.all1.gl25))][33:48]),
#        code = 3, angle = 90, length = 0.05, col = 2)
legend("bottomleft", expression(bold(B)), bty = "n", cex = 1.3)
abline(h = rli.all1.gl25[2, 1], lty = 2)
legend("bottomright", c("Group-specific", "Fixed"),
       lty = c(3, 0), pch=c(NA, 19),
       bty = "n", lwd=2)
dev.off()

rm(list=ls())
###

