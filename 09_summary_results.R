#### SUMMARIZING THE RESULTS ####
all.crit <- readRDS('data/all.criteria.rds')
all.crit2 <- readRDS('data/all.criteria_herb.rds')
sc_trees <- read.csv('sc_trees_list.csv', sep = ';', colClasses = c('NULL', NA, NA))

library(webr)
library(ggplot2)
library(cowplot)
library(multcompView)
library(circlize)
library(tidyverse)
library(ggpubr)

##### THE CONSERVATION STATUS OF SANTA CATARINA FLORA #####


#How many species were assessed?
dim(all.crit)[1] # 502 species in inventories only
dim(all.crit2)[1] # 555 species in inventories+herbaria

#How many threatened species?
all.crit <- all.crit %>% dplyr::rename(category.inv = category)
all.crit2 <- all.crit2 %>% dplyr::rename(category.herb = category)
names(all.crit2)[39] <- "category.inv"

## Overall
sum(table(all.crit$category.inv)[c("CR", "EN", "VU")]) # 443 threatened INVENTORY
sum(table(all.crit2$category.herb)[c("CR", "EN", "VU")]) # 477 threatened HERBARIA

#Threatened %
round(100 * sum(table(all.crit$category.inv)[c("CR", "EN", "VU")]) / dim(all.crit)[1], 1) #88.2% INVENTORY
round(100 * sum(table(all.crit2$category.herb)[c("CR", "EN", "VU")]) / dim(all.crit2)[1], 1) #85.9% HERBARIA


#Red List Index for Santa Catarina
red::rli(all.crit$category.inv, boot = TRUE, runs = 50000) #SC 0.5139 [0.4988-0.5294] INVENTORY
red::rli(all.crit2$category.herb, boot = TRUE, runs = 50000) #SC 0.534 [0.52-0.549] HERBARIA

#More frequent categories and criteria
round(100 * table(all.crit$category.inv)[c("CR", "EN", "VU", "NT", "LC")] / 
        dim(all.crit)[1], 1) #11.4% CR, 33.5% EN, 43.4% VU INVENTORY

round(100 * table(all.crit2$category.herb)[c("CR", "EN", "VU", "NT", "LC")] / 
        dim(all.crit)[1], 1) #8.2% CR, 39.4% EN, 47.4% VU HERBARIA


#How many near threatened
round(100 * sum(table(all.crit$category.inv)[c("NT")]) / dim(all.crit)[1], 1) #10.4%

round(100 * sum(table(all.crit2$category.herb)[c("NT")]) / dim(all.crit)[1], 1) #11.4%

#IUCN categories INVENTORY
tmp0 <- table(all.crit$main.criteria[!all.crit$category.inv %in% c("LC","NA","DD","NT")])
tmp <- sort(round(100 *  tmp0 / 
                    dim(all.crit[!all.crit$category.inv %in% c("LC","NA","DD","NT"),])[1], 4))
sum(tmp[grepl("A1", names(tmp))]) #0% of the threatened species had A1 within its main criteria INVENTORY
sum(tmp[grepl("A2", names(tmp))]) #77.87% of the threatened species had A2 within its main criteria INVENTORY
sum(tmp[grepl("B2", names(tmp))]) #34.53% of the threatened species had B2 within its main criteria INVENTORY
sum(tmp[grepl("B1", names(tmp))]) #12.86% of the threatened species had B1 within its main criteria INVENTORY
sum(tmp[grepl("C1", names(tmp))]) #No threatened species had C within its main criteria INVENTORY
sum(tmp[grepl("C2", names(tmp))]) #No threatened species had C within its main criteria INVENTORY
sum(tmp[grepl("D", names(tmp))]) #No threatened species had D within its main criteria INVENTORY

#IUCN categories HERBARIA
tmp0 <- table(all.crit2$main.criteria[!all.crit2$category.herb %in% c("LC","NA","DD","NT")])
tmp <- sort(round(100 *  tmp0 / 
                    dim(all.crit2[!all.crit2$category.herb %in% c("LC","NA","DD","NT"),])[1], 4))
sum(tmp[grepl("A1", names(tmp))]) #0% of the threatened species had A1 within its main criteria HERBARIA
sum(tmp[grepl("A2", names(tmp))]) #80.71% of the threatened species had A2 within its main criteria HERBARIA
sum(tmp[grepl("B2", names(tmp))]) #29.55% of the threatened species had B2 within its main criteria HERBARIA
sum(tmp[grepl("B1", names(tmp))]) #12.57% of the threatened species had B1 within its main criteria HERBARIA
sum(tmp[grepl("C1", names(tmp))]) #0.20% threatened species had C within its main criteria HERBARIA
sum(tmp[grepl("C2", names(tmp))]) #No threatened species had C within its main criteria HERBARIA
sum(tmp[grepl("D", names(tmp))]) #No threatened species had D within its main criteria HERBARIA


#Mean population reduction (won't change from inv to herbs — can only be measured for inventories)
summary(all.crit$reduction_A12); hist(all.crit$reduction_A12, nclass=40) #51.50% mean reduction
table(all.crit$reduction_A12>=90) #No species
table(all.crit$reduction_A12>=80) #No species
table(all.crit$reduction_A12>=70) #12 species, 2%
table(all.crit$reduction_A12>=50) #343, 68%
table(all.crit$reduction_A12<=30) #31, 6%

#Last 25 years
table(all.crit$reduction_A12.25ys>=90) #No species
table(all.crit$reduction_A12.25ys>=80) #No species
table(all.crit$reduction_A12.25ys>=70) #No species
table(all.crit$reduction_A12.25ys>=50) #2 species,0.3%
table(all.crit$reduction_A12.25ys<=30) #111 species, 22%


##### EOO and AOO ##### INVENTORY
#Filtering only threatened categories
tmp = all.crit %>% filter(category.inv %in% c("CR", "VU", "EN"))

#SMALL AOO AND EOO
summary(all.crit$AOO)
hist(log10(all.crit$AOO), nclass=40)
abline(v=log10(c(10,500,2000)), col=c("red","darkorange","gold"))
confint(lm(all.crit$AOO ~ 1))
summary(all.crit$nbe_loc_total)
hist(log10(all.crit$nbe_loc_total), nclass=40)
abline(v=log10(c(1,5,10)), col=c("red","darkorange","gold"))

confint(lm(all.crit$nbe_loc_total ~ 1))
100*table(all.crit$nbe_loc_total>10)/dim(all.crit)[1]
100*table(all.crit$sever.frag)/dim(all.crit)[1]

##### EOO and AOO ##### HERBARIA
#AOO
summary(all.crit2$AOO)# HERBARIA
hist(log10(all.crit2$AOO), nclass=40)# HERBARIA
abline(v=log10(c(10, 500, 2000)), col=c("red","darkorange","gold"))# HERBARIA
confint(lm(all.crit2$AOO ~ 1))# HERBARIA
table(all.crit2$AOO < 2000 & all.crit2$nbe_loc_total <= 10) #Threatened under AOO (28.28%)# HERBARIA
table(all.crit2$AOO < 2000 & all.crit2$nbe_loc_total <=10 & # HERBARIA
        all.crit2$sever.frag == TRUE) #Threatened under AOO and severely fragmented (19.63%)# HERBARIA

#EOO
summary(all.crit2$EOO)# HERBARIA
hist(log10(all.crit2$EOO), nclass=40)# HERBARIA
abline(v=log10(c(10, 500, 2000)), col=c("red","darkorange","gold"))# HERBARIA
confint(lm(all.crit2$EOO ~ 1))# HERBARIA
table(all.crit2$EOO < 20000 & all.crit2$nbe_loc_total <= 10) #Threatened under EOO (17.83%)# HERBARIA
table(all.crit2$EOO < 20000 & all.crit2$nbe_loc_total &# HERBARIA
        
        all.crit2$sever.frag == TRUE) #Threatened under AOO and severely fragmented (11.71%)# HERBARIA


#####################
#NUMBER OF LOCATIONS#
#####################
#How many species occurred in more than 10 locations? INVENTORY
summary(all.crit$nbe_loc_total) #INVENTORY
hist(log10(all.crit$nbe_loc_total), nclass=40) #INVENTORY
abline(v=log10(c(1, 5, 10)), col=c("red","darkorange","gold")) #INVENTORY
confint(lm(all.crit$nbe_loc_total ~ 1)) #INVENTORY
table(all.crit$nbe_loc_total>10)/502 #66.73% occurs in more than 10 localities #INVENTORY
table(all.crit$sever.frag>10)/502 #No species that occur in more than 10 localities are severely fragmented

#How many species are severely fragmented? #INVENTORY
table(all.crit$sever.frag) #25.09% of all species are severely fragmented #INVENTORY


#####################
#NUMBER OF LOCATIONS#
#####################
#How many species occurred in more than 10 locations? HERBARIA
summary(all.crit2$nbe_loc_total) #HERBARIA
hist(log10(all.crit2$nbe_loc_total), nclass=40) #HERBARIA
abline(v=log10(c(1, 5, 10)), col=c("red","darkorange","gold")) #HERBARIA
confint(lm(all.crit2$nbe_loc_total ~ 1)) #HERBARIA
table(all.crit2$nbe_loc_total>10)/555 #71.71% occurs in more than 10 localities
table(all.crit2$sever.frag>10)/555 #No species that occur in more than 10 localities are severely fragmented

#How many species are severely fragmented? #HERBARIA
table(all.crit2$sever.frag) #20.72% of all species are severely fragmented #HERBARIA



##### Main changes between categories and type of data (inv x inv+herbaria) #####
df = all.crit2[, c("category.herb", "category.inv")]
df[is.na(df)] <- "DD"

#Now we need to bind missing species from the assessment
#From 577 species, only 555 were assessed on inv+herb and 502 on inv
#Thus, we need to add 75 DD species to inv and 22 NE to herb
missing_sp <- data.frame(category.herb = rep("NE", 22), category.inv = rep("DD", 22))
df <- rbind(df, missing_sp)


#Inventory vs Herbaria
tmp <- table(paste0(df$category.herb,"_herb"), df$category.inv)
#LC to CR
100 * sum(tmp[c("CR_herb"), c("LC","NT")]) / sum(tmp) #no species became more threatened


#Threatened to not Threatened
100 * sum(tmp[c("LC_herb", "NT_herb"), c("CR", "EN", "VU")]) /
  sum(tmp) #1.9% of species became not threatened

count(df %>% filter(category.herb != "NE" & category.inv == "DD")) #53 species were assessed
#after adding herbaria data

# Occurrence insided protected areas
table(all.crit$protected>0) / dim(all.crit)[1] #65.33% of occurrences are inside PAs INVENTORY
table(all.crit2$protected>0) / dim(all.crit2)[1] #70.99% of occurrences are inside PAs HERBARIA

#INVENTORY
table(all.crit$protected[all.crit$category.inv %in% c("CR","EN","VU")]>0) / 
  dim(all.crit[all.crit$category.inv %in% c("CR","EN","VU"),])[1] #62.97% of threat. occurr. are inside PAs

#HERBARIA
table(all.crit2$protected[all.crit2$category.herb %in% c("CR","EN","VU")]>0) / 
  dim(all.crit2[all.crit2$category.herb %in% c("CR","EN","VU"),])[1] #68.13% of threat. occurr. are inside PAs


#INVENTORY
#Analyses of criticall endangered (CR) only
cr_spp = all.crit %>% filter(category.inv == "CR")
table(cr_spp$main.criteria) #2 under EOO, 55 under AOO
mean(cr_spp$reduction_A12) ##Mean reduction of 49.90% in the last 25 years
mean(cr_spp$reduction_A12.25ys) #Mean reduction of 33.17% in the last 25 years
table(cr_spp$protected) #94% of all CR species are not inside PAs

#HERBARIA
#Analyses of criticall endangered (CR) only
cr_spp = all.crit2 %>% filter(category.herb == "CR")
table(cr_spp$main.criteria) #2 under EOO, 35 under AOO, 4 under AOO+EOO
table(cr_spp$protected) #85.36% of all CR species are not inside PAs



##### THE INFLUENCE OF USING DIFFERENT CRITERIA ##### INVENTORY #####
##Species classified under all criteria vs. criterion A and criterion B
tmp <- all.crit[!is.na(all.crit$reduction_A12) &  # Assessed under criterion A
                  !is.na(all.crit$AOO) & # Assessed under criterion B2
                  !all.crit$category.inv %in% "NA", ]
dim(tmp)[1] ## 502 species
100*dim(tmp)[1]/dim(all.crit)[1] #100% of species classified within criteria A or B
tmp$category.inv[tmp$category.inv %in% c("LC","NT")] <- "LC+NT"

#All species
mat1 <- table(tmp[, c("category.inv","category_A")]) ## All vs. A
mat2 <- table(tmp[, c("category.inv","category_B")]) ## All vs. B
mat3 <- table(tmp[, c("category.inv","category_C")]) ## All vs. C
mat4 <- table(tmp[, c("category.inv","category_D")]) ## All vs. D
mat5 <- table(tmp[, c("category_A","category_B")]) ## A vs. B

#Red List Index per category
rl0 <- round(red::rli(tmp$category.inv, boot = TRUE, runs = 50000), 4) #0.46 0.473 0.486
rl1 <- round(red::rli(tmp$category_A, boot = TRUE, runs = 50000), 4) #0.446 0.454 0.462
rl2 <- round(red::rli(tmp$category_B, boot = TRUE, runs = 50000), 4) #0.346 0.367 0.389
rl3 <- round(red::rli(tmp$category_C, boot = TRUE, runs = 10),4) #0.567 0.581 0.6
rl4 <- round(red::rli(tmp$category_D, boot = TRUE, runs = 50000), 4) #0.6 0.6 0.6

#% of threatened species per criteria
ts0 <- round(100*sum(mat1[-3,]) / sum(mat1),2) #88.25%, all criteria
ts1 <- round(100*sum(mat1[,-3]) / sum(mat1),2) #74.5%, only A
ts2 <- round(100*sum(mat2[,-3]) / sum(mat2),2) #31.67%, only B
ts3 <- round(100*sum(mat3[,-2]) / sum(mat3),2) #1.59%, only C
ts4 <- round(100*sum(mat4[,-1]) / sum(mat4),2) #16.53%, only D


## Table (S)1? ##
tabS1 <- matrix(c(ts1, paste0(rl1[2]," [",paste0(round(rl1[c(1,3)],3), collapse = "-"),"]"),
                  ts2, paste0(rl2[2]," [",paste0(round(rl2[c(1,3)],3), collapse = "-"),"]"),
                  ts3, paste0(rl3[2]," [",paste0(round(rl3[c(1,3)],3), collapse = "-"),"]"),
                  ts4, paste0(rl4[2]," [",paste0(round(rl4[c(1,3)],3), collapse = "-"),"]"),
                  ts0, paste0(rl0[2]," [",paste0(round(rl0[c(1,3)],3), collapse = "-"),"]")),
                  nrow = 5, ncol = 2, byrow = TRUE,
                  dimnames = list(c("catA","catB","catC","catD","all"), c("Threatened","RLI")))
tabS1
write.csv(tabS1,"TableS1.csv")



##### THE INFLUENCE OF USING DIFFERENT CRITERIA ##### HERBARIA #####
##Species classified under all criteria vs. criterion A and criterion B
tmp <- all.crit2[!is.na(all.crit2$reduction_A12) &  # Assessed under criterion A
                  !is.na(all.crit2$AOO) & # Assessed under criterion B2
                  !all.crit2$category.herb %in% "NA", ]
dim(tmp)[1] ## 502 species
100*dim(tmp)[1]/dim(all.crit2)[1] #90.45% of species classified within criteria A or B
tmp$category.herb[tmp$category.herb %in% c("LC","NT")] <- "LC+NT"

#All species
mat1 <- table(tmp[, c("category.herb","category_A")]) ## All vs. A
mat2 <- table(tmp[, c("category.herb","category_B")]) ## All vs. B
mat3 <- table(tmp[, c("category.herb","category_C")]) ## All vs. C
mat4 <- table(tmp[, c("category.herb","category_D")]) ## All vs. D
mat5 <- table(tmp[, c("category_A","category_B")]) ## A vs. B

#Red List Index per category
rl0 <- round(red::rli(tmp$category.inv, boot = TRUE, runs = 50000), 4) #0.498 0.514 0.529
rl1 <- round(red::rli(tmp$category_A, boot = TRUE, runs = 50000), 4) #0.446 0.454 0.462
rl2 <- round(red::rli(tmp$category_B, boot = TRUE, runs = 50000), 4) #0.396 0.418 0.441
rl3 <- round(red::rli(tmp$category_C, boot = TRUE, runs = 10),4) #0.529 0.563 0.6
rl4 <- round(red::rli(tmp$category_D, boot = TRUE, runs = 50000), 4) #0.6 0.6 0.6

#% of threatened species per criteria
ts0 <- round(100*sum(mat1[-3,]) / sum(mat1),2) #86.06%, all criteria #ajeitar na mão WTF
ts1 <- round(100*sum(mat1[,-3]) / sum(mat1),2) #74.5%, only A
ts2 <- round(100*sum(mat2[,-3]) / sum(mat2),2) #21.71%, only B
ts3 <- round(100*sum(mat3[,-2]) / sum(mat3),2) #1.59%, only C
ts4 <- round(100*sum(mat4[,-1]) / sum(mat4),2) #16.53%, only D


## Table (S)1? ##
tabS1 <- matrix(c(ts1, paste0(rl1[2]," [",paste0(round(rl1[c(1,3)],3), collapse = "-"),"]"),
                  ts2, paste0(rl2[2]," [",paste0(round(rl2[c(1,3)],3), collapse = "-"),"]"),
                  ts3, paste0(rl3[2]," [",paste0(round(rl3[c(1,3)],3), collapse = "-"),"]"),
                  ts4, paste0(rl4[2]," [",paste0(round(rl4[c(1,3)],3), collapse = "-"),"]"),
                  ts0, paste0(rl0[2]," [",paste0(round(rl0[c(1,3)],3), collapse = "-"),"]")),
                nrow = 5, ncol = 2, byrow = TRUE,
                dimnames = list(c("catA","catB","catC","catD","all"), c("Threatened","RLI")))
tabS1
write.csv(tabS1,"TableS1_herb.csv")


###################
##### FIGURES #####
###################
source('scripts/myPieDonut.R')

{
## Colors for each category
cores <- c(CR = "red", EN = "darkorange", VU = "gold", 
           NT = "yellowgreen", LC = "forestgreen", `NA` = "lightgrey", NE = "darkgrey")


##################
#### FIGURE 1 ####
##################
#Creating and organizing the data frame to be plotted
pie.df <- all.crit[,c("category.inv","main.criteria")]
names(pie.df)[1] <- "category"
pie.df$category[pie.df$category %in% "CR_PE"] <- "CR" 
pie.df$main.criteria[pie.df$main.criteria %in% "A1+A2+B1+B2+C1+C2+D"] <- "all" 
pie.df$main.criteria[pie.df$main.criteria %in% c("A1+A2+B1+B2+C1","A1+A2+B2+C1")] <- "other" #"A+B+C"
pie.df$main.criteria[pie.df$main.criteria %in% c("A2+C1","A1+A2+C1+C2","A1+A2+C1")] <- "other" #"A+C"
pie.df$main.criteria[pie.df$main.criteria %in% c("A1+A2+B2","A2+B1+B2","A2+B2","A1+A2+B1+B2")] <- "A2, B2"
pie.df$main.criteria[pie.df$main.criteria %in% c("B1")] <- "other"
pie.df$main.criteria[pie.df$main.criteria %in% c("A1+A2+C2+D")] <- "other" #"A+C+D"
pie.df$main.criteria[is.na(pie.df$main.criteria)] <- "other"
pie.df$main.criteria[pie.df$main.criteria %in% c("A1+A2")] <- "A2"
pie.df$main.criteria[pie.df$category %in% "NA"] <- ""
pie.df <- pie.df[order(match(pie.df$category, names(cores)[-1])),]
pie.df <- pie.df[!pie.df$category %in% "DD",]
pie.df$main.criteria[pie.df$main.criteria %in% c("B1+B2")] <- "B1+2"

#### PIE DONUT #### INVENTORY
dados <- pie.df

pie.all <- my.PieDonut(dados, aes(category, main.criteria),
                       ratioByGroup=FALSE, showPieName = FALSE,
                       start=75,r0=0, r1=1,r2=1.3,
                       showRatioThreshold = 0.02, labelpositionThreshold =  0.03,
                       main.colors = cores,
                       #donut.colors = NULL,
                       pieAlpha = 0.75,
                       plot.piedonut = FALSE,
                       return = "pie", pieLabelSize = 4, donutLabelSize = 3,
                       title = "A) Forest surveys only")

donut.all <- my.PieDonut(dados, aes(category, main.criteria),
                         ratioByGroup=FALSE, showPieName = FALSE,
                         start=75,r0=0, r1=1,r2=1.3,
                         showRatioThreshold = 0.02, labelpositionThreshold =  0.03,
                         main.colors = cores,
                         #donut.colors = NULL,
                         pieAlpha = 0.75,
                         plot.piedonut = FALSE,
                         return = "donut", pieLabelSize = 4, donutLabelSize = 3,
                         title = "A) Forest surveys only")

## Constructing the figures
pie.donut.all <- ggdraw(pie.all) + draw_plot(donut.all)
pie.donut.all  

#Creating and organizing the data frame to be plotted #### HERBARIA
pie.df <- all.crit2[,c("category.herb","main.criteria")]
names(pie.df)[1] <- "category"
pie.df$category[pie.df$category %in% "CR_PE"] <- "CR" 
pie.df$main.criteria[pie.df$main.criteria %in% "A1+A2+B1+B2+C1+C2+D"] <- "all" 
pie.df$main.criteria[pie.df$main.criteria %in% c("A1+A2+B1+B2+C1","A1+A2+B2+C1")] <- "other" #"A+B+C"
pie.df$main.criteria[pie.df$main.criteria %in% c("A2+C1","A1+A2+C1+C2","A1+A2+C1")] <- "other" #"A+C"
pie.df$main.criteria[pie.df$main.criteria %in% c("A1+A2+B2","A2+B1+B2","A2+B2","A1+A2+B1+B2")] <- "A2, B2"
pie.df$main.criteria[pie.df$main.criteria %in% c("B1")] <- "other"
pie.df$main.criteria[pie.df$main.criteria %in% c("A1+A2+C2+D")] <- "other" #"A+C+D"
pie.df$main.criteria[is.na(pie.df$main.criteria)] <- "other"
pie.df$main.criteria[pie.df$main.criteria %in% c("A1+A2")] <- "A2"
pie.df$main.criteria[pie.df$category %in% "NA"] <- ""
pie.df <- pie.df[order(match(pie.df$category, names(cores)[-1])),]
pie.df <- pie.df[!pie.df$category %in% "DD",]
pie.df$main.criteria[pie.df$main.criteria %in% c("B1+B2")] <- "B1+2"

#### PIE DONUT #### HERBARIA
dados <- pie.df

pie.all <- my.PieDonut(dados, aes(category, main.criteria),
                       ratioByGroup=FALSE, showPieName = FALSE,
                       start=75,r0=0, r1=1,r2=1.3,
                       showRatioThreshold = 0.02, labelpositionThreshold = 0.03,
                       main.colors = cores,
                       #donut.colors = NULL,
                       pieAlpha = 0.75,
                       plot.piedonut = FALSE,
                       return = "pie", pieLabelSize = 4, donutLabelSize = 3,
                       title = "B) Forest surveys plus herbaria")

donut.all <- my.PieDonut(dados, aes(category, main.criteria),
                         ratioByGroup=FALSE, showPieName = FALSE,
                         start=75,r0=0, r1=1,r2=1.3,
                         showRatioThreshold = 0.02, labelpositionThreshold =  0.03,
                         main.colors = cores,
                         #donut.colors = NULL,
                         pieAlpha = 0.75,
                         plot.piedonut = FALSE,
                         return = "donut", pieLabelSize = 4, donutLabelSize = 3,
                         title = "B) Forest surveys plus herbaria")

## Constructing the figures
pie.donut.all2 <- ggdraw(pie.all) + draw_plot(donut.all)
pie.donut.all2 


final = ggarrange(pie.donut.all, pie.donut.all2,
                  ncol = 2,
                  font.label = list(size=16),
                  align = "h")

final = annotate_figure(final, 
                        top = text_grob("Comparison between IUCN categories",
                        face = "bold", size = 16))

ggsave2("Figure1_NEW_herb.jpg", final, "jpeg", "figures/",
        width = 30, height = 20, units = "cm", dpi = 300)
}
##### FIGURE 2 #####

##### NEW DATA FOR IMAGES #####
#Considering all possible species that may occur in SC

#Species that are going to be considered data deficient
#Remove subsp...
del <- c("Inga vera subsp. affinis", "Piptocarpha axillaris subsp. axillaris")

not_found_inv <- setdiff(sc_trees$species, all.crit$species) #Not Evaluated
not_found_inv <- not_found_inv[!not_found_inv %in% del]
not_found_herb = setdiff(sc_trees$species, all.crit2$species) #Not Evaluated
not_found_herb <- not_found_herb[!not_found_herb %in% del]

not_found_inv <- read.csv('C:/Users/Public/Desktop/not_found_inv.csv', sep = ';')
not_found_herb <- read.csv('C:/Users/Public/Desktop/not_found_herb.csv', sep = ';')

all.crit <- rbind(all.crit, not_found_inv)
all.crit2 <- rbind(all.crit2, not_found_herb)



# Previous vs. new assessments
library(circlize)

####Figure parameters####
{
jpeg(filename = "figures/Figure2_NEW_TEST.jpg", width = 3750, height = 2000, units = "px", pointsize = 12,
     res = 300, bg="white", quality = 100)

#### NEW VS. PREVIOUS (NATIONAL) ####
mat <- as.matrix(table(paste0(all.crit2$category.herb, "_new"), 
                       paste0(all.crit2$category.inv, "_prev")))
mat
mat <- mat[c(1,2,6,5,3,4), rev(c(1,2,6,5,3,4))]
mat

#Defining the colors of tracks and links
grid.col = c(CR_prev = "red", 
             EN_prev = "darkorange", 
             VU_prev = "gold", 
             NT_prev = "yellowgreen", 
             LC_prev = "forestgreen", 
             NE_prev = "lightgrey",
             CR_new = "red", 
             EN_new = "darkorange", 
             VU_new = "gold",
             NT_new = "yellowgreen",
             LC_new = "forestgreen",
             NE_new = "lightgrey"
)

col_mat = rep(rev(c("red", "darkorange", "gold", "yellowgreen", "forestgreen", "lightgrey")),
              each=6)
#col_mat[mat >= 15] = adjustcolor(col_mat[mat >= 15], alpha.f = 0.5)
#col_mat[mat < 15] = adjustcolor(col_mat[mat < 15], alpha.f = 0.9)
#col_mat[mat < 5] = "#00000000"

mat[mat < 10 & mat >= 5] = mat[mat < 10 & mat >= 5] + 1
mat[mat > 0 & mat < 5] = mat[mat > 0 & mat < 5] + 2
transp <- rep(0.3, length(col_mat))
transp[col_mat %in% "forestgreen"] <- 0.6


#plotting the diagram
circos.clear()
circos.par(start.degree = 90)
visible = matrix(TRUE, nrow = nrow(mat), ncol = ncol(mat))
#diag(visible) = FALSE
#lava::revdiag(visible[,-c(1:2)]) = FALSE
chordDiagram(mat, big.gap = 10, annotationTrack = "grid", annotationTrackHeight = mm_h(5),
             grid.col = grid.col, col = col_mat,
             transparency = transp,
             self.link = 1, #link.visible = visible,
             h=1.1, #w=0.5,
             #direction.type = "arrows", link.arr.length = 0.2, link.arr.width = 0.1, directional = -1,
             link.lwd = 4,
             h.ratio = 0.9,
             #reduce_to_mid_line = TRUE,
             w2=0.5, rou=0.2
             #point1 = rep(0,16)
             )

#Putting legends on
sec.ind <- c("CR","EN","VU","NT","LC","NE","NE","LC","NT","VU","EN","CR")#,"EW","EX")
for(si in get.all.sector.index()) {
  lab <- sec.ind[which(si == get.all.sector.index())]
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
    if(si == "VU_hi") {
      circos.text(mean(xlim), mean(ylim), labels = "VU", 
                  sector.index = si, track.index = 1, cex = 1.1, adj= 0.1,
                  facing = "inside", niceFacing = FALSE, col = "black")
      
    } else {
  circos.text(mean(xlim), mean(ylim), labels = lab, 
              sector.index = si, track.index = 1, cex = 1.1, #adj= 0.1,
              facing = "inside", niceFacing = TRUE, col = "black")
    }  
}
legend("topleft","Forest inventory assess.", bty="n", cex=1.2, adj=c(-0.55,3))
legend("topright","Forest inventory plus herbaria assess.", bty="n", cex=1.2, adj=c(0.10,3))
legend("top",legend=expression(bold("Changes in IUCN categories after adding herbaria data")),
       bty="n",horiz=F,cex=1.5,x.intersp=-0.7,y.intersp=-0.3)

dev.off()
}
  
