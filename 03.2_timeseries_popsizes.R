#Load data
fit_lower <- readRDS('data/02_lower_popsize.rds')
fit_upper <- readRDS('data/02_fit_popsize.rds')
fit_fit <- readRDS('data/02_upper_popsize.rds')
fit_fit <- as.data.frame(fit_fit)

#Loading trees per cell table
trees_lower <- readRDS('data/01_Ntrees_lower.rds')
trees_fit <- readRDS('data/01_Ntrees_fit.rds')
trees_upper <- readRDS('data/01_Ntrees_upper.rds')

#Separate a df for each information
#lower bound of population size
spp_lower <- fit_lower[, 4:505]
years_lower <- trees_lower[, 8:43]
#fit bound of population size
spp_fit <- fit_fit[, 4:505]
years_fit <- trees_fit[, 8:43]
#upper bound of population size
spp_upper <- fit_upper[, 4:505]
years_upper <- trees_upper[, 8:43]

#Creating a matrix output based on the dimensions of number of species (502) and number of years (35)
#Do the same for the upper and lower CI
output <- matrix(0, dim(spp_upper)[2], dim(years_upper)[2])
colnames(output) <- names(years_upper)
rownames(output) <- names(spp_upper)

#Looping: summing the values for each year and saving
for (y in names(years_upper)){
  for (sp in names(spp_upper)){
    Sum_nTreebySp <- sum(years_upper[, y] * spp_upper[, sp])
    output[sp, y] <- Sum_nTreebySp 
  }
}

#Dataframe
output_df_upper <- as.data.frame(output) 
sum(output_df_upper$trees_2019) #~90% of all trees

#Saving
write.csv(output_df_upper, 'timeseries_upper_popsize.csv')
saveRDS(output_df_upper, 'data/03_timeseries_upper_popsize.rds')


rm(list=ls())

