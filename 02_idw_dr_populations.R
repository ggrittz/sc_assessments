library(raster)
library(gstat)
library(data.table)
library(tidyverse)

####Start####
rm(list = ls())

#Loading RDS
points = readRDS('data/00_points_filtered_nm.rds')
grid_env = readRDS('data_old_stuff/00_grid_env.rds')
sc_wo_grassland = readRDS('data_old_stuff/00_sc_wo_grassland.rds')

####Data adjustments####
#Species list from abundance_data_set_per_sitecode_v1.csv
species_list = readRDS('data/checked_specieslist.rds')

#Filtering what is needed
species_list = species_list %>% filter(taxon.rank == "species")
sum(is.na(species_list))

#Subsetting coordinates and an ID from TreeCo points
plots_coord <- points %>% dplyr::select(SiteCode.x, long1, lat1)
plots_coord$SiteCode.x <- as.character(plots_coord$SiteCode.x)

#Obtaining coordinates point for each register based on the ID
joined_data <- merge(species_list, plots_coord, by.x = 'SiteCode', by.y = 'SiteCode.x', sort = FALSE, all.y = TRUE)
sum(is.na(joined_data))

#Choosing only what is a tree as defined before
trees = read.csv('sc_trees_list.csv', header = TRUE, sep = ';')
trees = trees$species

joined_data <- joined_data[joined_data$species.correct %in% trees]

#How many species?
spp = unique(joined_data$species.correct)

#Calculating correct N
joined_data$Ntotal_new = with(joined_data, ave(N, SiteCode, FUN = sum))
#Adjusting to effort
joined_data$N = joined_data$N / joined_data$effort_ha
joined_data$Ntotal_new = joined_data$Ntotal_new / joined_data$effort_ha

#Calculating relative abundance
setDT(joined_data)[, rel_ab:=N/Ntotal_new, by = SiteCode]

#Testing if sum = 1
teste = subset(joined_data, SiteCode %in% "SCcelso2")
sum(teste$rel_ab)

#Assign equal CRS to all spatial data
new_crs <- "+proj=utm +zone=22 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
coordinates(joined_data)=~long1+lat1
crs(joined_data)=crs(grid_env)=new_crs

#Species vector
species <- sort(unique(joined_data$species.correct))

#Creating an empty list to save the IDW results
sp_list <- setNames(vector('list', length(species)), species)

####Inverse Distance Weightning####
for (sp in species) {
  gstat_spec <- joined_data[joined_data$species.correct == sp, ]
  out <- idw(rel_ab ~ 1, gstat_spec, grid_env, nmax=50, maxdist=0.125, idp=2, debug.level=0)
  out@data[out@data < 0] <- 0
  out <- crop(out, sc_wo_grassland)
  sp_list[[sp]] <- out[[1]]
  print(sp)
}

#Binding the results
sp_list2<-do.call(cbind, sp_list)
sp_list3 <- as.data.frame(sp_list2)
sp_list3[is.na(sp_list3)] = 0
sum_populations <- sum(sp_list3, na.rm=T)
sum_populations

#Joining coordinates/trees per ha for each fit of the kriging model (lower, fit, upper)
vegtable <- read.csv('sc_wo_grassland_df.csv', header = TRUE, sep = ',')

df_fit <- cbind(vegtable$X, vegtable$Y, vegtable$krig_fit, sp_list3)
names(df_fit)[1] <- "longitude"
names(df_fit)[2] <- "latitude"
names(df_fit)[3] <- "krig_fit"

df_lower <- cbind(vegtable$X, vegtable$Y, vegtable$krig_lower, sp_list3)
names(df_lower)[1] <- "longitude"
names(df_lower)[2] <- "latitude"
names(df_lower)[3] <- "krig_lower"

df_upper <- cbind(vegtable$X, vegtable$Y, vegtable$krig_upper, sp_list3)
names(df_upper)[1] <- "longitude"
names(df_upper)[2] <- "latitude"
names(df_upper)[3] <- "krig_upper"

#Testing a plot
coordinates(df_upper)<-~longitude+latitude
spplot(df_upper, 'Araucaria.angustifolia')

#Saving what is needed
write.csv(as.data.frame(joined_data), 'spp_coords.csv')
#Saving RDS
saveRDS(sp_list3, 'data/02_specieslist_DR.rds')
saveRDS(df_lower, 'data/02_lower_popsize.rds')
saveRDS(df_fit, 'data/02_fit_popsize.rds')
saveRDS(df_upper, 'data/02_upper_popsize.rds')
####End####
rm(list = ls())
