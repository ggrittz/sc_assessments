#Start
rm(list = ls())

#####Libraries####
{
  library(sf)
  library(raster)
  library(rgdal)
  library(gstat)
  library(sp)
  library(dplyr)
  library(automap)
  library(rgeos)
  library(ggplot2)
  library(tidyverse)
}

#####Loading and adjusting data####
sc_buffer = readRDS('rds_files/00_sc_buffer.rds')
#saveRDS(sc_buffer, 'data/00_sc_buffer.rds')
sc_wo_grassland = readRDS('rds_files/00_sc_wo_grassland.rds')
#saveRDS(sc_wo_grassland, 'data/00_sc_wo_grassland.rds')
grid_env = readRDS('rds_files/00_grid_env.rds')
#saveRDS(grid_env, 'data/00_grid_env.rds')

#Southern Brazil TreeCo points
points = read.csv('treeco/treeco_for_guilherme_grittz/sites_data_set_v1.csv',
                  header = T, sep = ',')#Cropping all points that are inside the state's buffer
coordinates(points) <- ~long1+lat1

#Points that are in the buffer area only
points = crop(points, sc_buffer)
plot(sc_buffer)
plot(points, add = TRUE)
dim(points)[1] #584 sample points inside Santa Catarina buffer

#Adjusting overall crs
crs(points) = crs(sc_wo_grassland) = crs(grid_env) = crs(sc_buffer)
points = as.data.frame(points)

#Any non-unique plots?
n_occur <- data.frame(table(points$SiteCode.x))
n_occur[n_occur$Freq > 1,]
dup = points[points$SiteCode.x %in% n_occur$Var1[n_occur$Freq > 1],]


#All species found in 'points'
species_list = read.csv('treeco/treeco_for_guilherme_grittz/abundance_data_set_per_sitecode_v1.csv',
                        header = T, sep = ',')

#Removing some NA columns
species_list <- species_list[,-which(names(species_list) == "SubSite")]
species_list <- species_list[,-which(names(species_list) == "DoR")]
species_list <- species_list[,-which(names(species_list) == "family")]
species_list <- species_list[,-which(names(species_list) == "genus")]
sum(is.na(species_list)) #all good

#Obtaining coordinates point for each register based on the SiteCode
joined_data <- merge(species_list, points[, c('SiteCode.x', 'long1', 'lat1')], 
                     by.x = 'SiteCode', by.y = 'SiteCode.x', sort = FALSE, all.y = TRUE)
sum(is.na(joined_data))
unique(joined_data$SiteCode) #581 sample points with coordinates
unique(points$SiteCode.x)

#TreeCo points: filtering points based on their confiability
points = points %>% filter(confiabilidade %in% c('boa', 'exata', 'precisa')) #3 points removed

#Remove problematic (old_survey and no_status)
points = points %>% filter(problems == "ok")

#Spatial before saving
coordinates(points)=~long1+lat1
saveRDS(points, 'data/00_points_filtered.rds')

#End
rm(list = ls())
