rm(list=ls())

df <- read.csv('sc_wo_grassland_df.csv', header = TRUE, sep = ',')

#Upper limit of 1850 being 2500 ha, that is, fully forested
df$ha_1850 <- as.numeric(2500)
df <- df[, c(1:10, 46, 11:45)]

#Tree quantity per cell from upper, fit, and lower confidence interval bounds
pop_lower <- df$krig_lower * df[, 11:46]
pop_fit <- df$krig_fit * df[, 11:46]
pop_upper <- df$krig_upper * df[, 11:46]

#Re-jonining the data
trees_lower <- cbind(df[, 1:7], pop_lower)
trees_fit <- cbind(df[, 1:7], pop_fit)
trees_upper <- cbind(df[, 1:7], pop_upper)

#Function to substitute ha for trees
colClean <- function(x){colnames(x) <- gsub("ha", "trees", colnames(x)); x} 

trees_lower <- colClean(trees_lower)
trees_fit <- colClean(trees_fit)
trees_upper <- colClean(trees_upper)

#These objects are then multiplied by DR from each species to obtain population sizes
saveRDS(trees_upper, 'data/01_Ntrees_upper.rds')
saveRDS(trees_fit, 'data/01_Ntrees_fit.rds')
saveRDS(trees_lower, 'data/01_Ntrees_lower.rds')

#Saving df again, now corrected for everything
saveRDS(df, 'data/00_ecoreg_krig_ha_df.rds')

rm(list=ls())
