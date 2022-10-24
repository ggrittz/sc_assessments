rm(list = ls())

#Libraries
library(ConR)
library(sf)
library(dplyr)
library(data.table)
library(rnaturalearthdata)
library(flora)
library(raster)
library(rgdal)

##### Criteria B: IUCN #####

#Loading inventory data
sc = readRDS('data_old_stuff/00_sc_buffer.rds')
occ = read.csv('spp_coords.csv', as.is = TRUE)

#Loading herbaria data
herb = read.csv('plantR_herbarium_filtered.csv')

#Obtaining families from Flora do Brasil (to add to inventory data)
sp.unique = unique(occ$species.correct)
families = get.taxa(sp.unique)
occ = merge(occ, families[, c("original.search", "family")], by.x = "species.correct", 
            by.y = "original.search", all.x = TRUE, sort = FALSE)

#Adjusting inventory data to ConR format
occ_data = occ %>% dplyr::select(lat1, long1, species.correct, family)
occ_data = occ_data %>% rename(ddlat = lat1, ddlon = long1, tax = species.correct)
occ_data$ddlat <- round(occ_data$ddlat, 2)
occ_data$ddlon <- round(occ_data$ddlon, 2)
head(occ_data)

#Adjusting herbarium data to ConR format
herb = herb %>% dplyr::select(decimalLatitude.new1, decimalLongitude.new1, suggestedName, family)
names(herb) <- c("ddlat", "ddlon", "tax", "family")

inv_herb = rbind(occ_data, herb)


##### EXTENT OF OCCURRENCE (EOO) ##### runned on ConR 1.4!
  EOO.hull <- ConR::EOO.computing(inv_herb[, 1:3],
                                  method = "convex.hull",
                                  method.less.than3 = "not comp",
                                  export_shp = TRUE,
                                  exclude.area = FALSE, country_map = NULL,
                                  write_shp = FALSE, 
                                  write_results=FALSE, file.name = "EOO.hull", 
                                  parallel = TRUE, NbeCores = 6) 

#Extracting the EOO from the object
EOO <- do.call(rbind.data.frame, EOO.hull[grepl("EOO", names(EOO.hull))])
sp_names <- unique(inv_herb$tax)
sp_names <- as.data.frame(sp_names)
sp_names <- sp_names[order(sp_names[,'sp_names']), ]
dimnames(EOO) <- list(sp_names, "EOO")
saveRDS(EOO, "data/04_EOO_numeric_herb.rds")
#EOO <- EOO.hull[[1]]
#EOO1 <- EOO.hull1[[1]]

#Extracting the geometry from the object
shps <- EOO.hull[grepl("spatial.polygon", names(EOO.hull))]
for (i in 1:length(shps))
  slot(slot(shps[[i]], "polygons")[[1]], "ID") <- sp_names[!is.na(EOO$EOO)][i]
shps <- do.call(rbind, shps)
shps_df <- SpatialPolygonsDataFrame(shps, data.frame(tax = names(shps), row.names = names(shps)))
#shps_df <-EOO.hull[[2]]
shps_df$tax <- as.character(shps_df$tax)

#Inspecting an example
sp <- "Araucaria angustifolia"
plot(shps_df[sp,])
points(occ_data[occ_data$tax %in% sp, 2:1], pch=19)
plot(sc, add=TRUE, border = "grey")

#Saving
shps_df_sf <- sf::st_as_sf(shps_df)
saveRDS(shps_df, "data/04_spp_convexhull_polygons_herb.rds")
saveRDS(shps_df_sf, "data/04_spp_convexhull_polygons_sf_uncropped_herb.rds")


##### SUBPOPULATIONS ##### New functions already on ConR 2.0
devtools::install_github("gdauby/ConR@devel")

radius <- subpop.radius(inv_herb[, c(1:3)], factor.div = 10, quant.max = 0.9)
sum(is.na(radius))
#81 spp NA! #50 NA with herbaria data

#Getting genus- and family-specific radius for missing estimates
tmp <- merge(radius, inv_herb[, c('family', 'tax')], by = 'tax', all.x = TRUE, sort=FALSE)
tmp <- tmp[order(tmp$tax), ]
tmp$genus <- as.character(sapply(strsplit(tmp$tax, " "), function(x) x[1]))
#Median values from genus and family
gen.radius <- aggregate(as.double(tmp$radius), list(tmp$genus), median, na.rm = TRUE)
fam.radius <- aggregate(as.double(tmp$radius), list(tmp$family), median, na.rm = TRUE)
#Merge family and genus radius
radius.new <- merge(tmp, gen.radius, by.x = 'genus', by.y = "Group.1", all.x = TRUE, sort=FALSE)
radius.new <- merge(radius.new, fam.radius, by.x = 'family', by.y = "Group.1", all.x = TRUE, sort=FALSE)
radius.new$radius.final <- as.double(radius.new$radius)
radius.new$radius.final[is.na(radius.new$radius.final)] <- radius.new$x.x[is.na(radius.new$radius.final)]
radius.new$radius.final[is.na(radius.new$radius.final)] <- radius.new$x.y[is.na(radius.new$radius.final)]
radius.new <- radius.new[order(radius.new$tax), ]
inv_herb <- inv_herb[order(inv_herb$tax), ]

#Checking IDs
table(radius$tax %in% inv_herb$tax)
#hist(radius.new$radius.final, nclass = 40)

#Dropping 'radius' column to avoid bad coding
radius.new <- radius.new[, -c(4)]
#Assign the median km as a maximum dispersal distance
median(radius.new$radius.final, na.rm=T) #28km
radius.new$radius.final[is.na(radius.new$radius.final)] <- 27.269 #if still NA, use the overall median

#Getting the subpopulations
radius.new <- radius.new[, c("tax", "radius.final")]
names(radius.new) <- c("tax", "radius")
saveRDS(radius.new, 'data/04_species_subpop_radius_herb.rds')
write.csv(radius.new, 'species_subpop_radius_herb.csv')

#Calculating subpopulations based on species, genus and family 1/10 max distance (90% quantile)
SUB1 <- subpop.comp(inv_herb[, c(1:3)],
                    Resol_sub_pop = radius.new[, c("tax", "radius")])

SUB1_teste <- as.data.frame(SUB1)
names(SUB1_teste)[1] <- "Number_subpop"

saveRDS(SUB1_teste, 'data/04_n_subpop_herb.rds')


##### AREA OF OCCUPANCY (AOO) #####

#Loading data
sc = readRDS('data_old_stuff/00_sc_buffer.rds')

#Occurrence data from above (occ_data)

#ConR 1.3 yet
AOO <- AOO.computing(inv_herb[, c(1:3)], 
                     Cell_size_AOO = 2, # Size of the grid cell in kilometers
                     nbe.rep.rast.AOO = 30, # number of raster with random starting position for estimating the AOO
                     parallel = TRUE, NbeCores = 6, 
                     show_progress = TRUE, export_shp = FALSE)

#Obtaining the df
AOO <- as.data.frame(AOO)
saveRDS(AOO, 'data/04_AOO_herb.rds')
gc()

                                 
##### NUMBER OF LOCATIONS ##### #devtools::install_github("gdauby/ConR@devel") NEW FUNCTIONS on ConR 2.0!
ucs <- readRDS('data_old_stuff/04_strict_ucs_sc.rds')

locations <- locations.comp(inv_herb[, c(1:3)], 
                            method = "fixed_grid",
                            nbe_rep = 30, # number of raster with random starting position
                            Cell_size_locations = 10, #grid size in kilometres used for estimating the number of location
                            Rel_cell_size = 0.05,
                            protec.areas = ucs, #SpatialPolygonsDataFrame, shapefile with protected areas.
                            ID_shape_PA= "NAME", #field name of protected areas with ID
                            method_protected_area="no_more_than_one", 
                            parallel = TRUE, NbeCores = 6)

table(names(locations[[3]]) == names(locations[[4]]))

tmp_loc <- cbind.data.frame(locations[[3]], locations[[4]])
names(tmp_loc) <- c("InsideUC", "OutsideUC")

tmp_loc <- tibble::rownames_to_column(tmp_loc, "species")
tmp_loc$species <- gsub("_", " ", tmp_loc$species)
tmp_loc <- tmp_loc[order(tmp_loc$species), ]
saveRDS(tmp_loc, 'data/04_number_of_locations_herb.rds')


##### SEVERE FRAGMENTATION #####
rm(list=ls())

#Getting the population radius for each species
tmp2 <- readRDS("data/04_species_subpop_radius_herb.rds")
tmp2$radius <- as.double(tmp2$radius)

#Getting the AOO for each species
AOO <- readRDS("data/04_AOO_herb.rds")
AOO <- tibble::rownames_to_column(AOO, "tax")

toto1 <- inv_herb[!is.na(inv_herb$ddlat) & !is.na(inv_herb$ddlon), c(1:3), ]
toto1 <- toto1[order(toto1$tax),]
table(toto1$tax == tmp2$tax)
toto1$radius <- tmp2$radius
toto1$radius <- round(toto1$radius, 0)
toto1 <- merge(toto1, AOO[, c("tax", "AOO")], by = "tax", all.x = TRUE, sort = FALSE)
toto1 <- toto1[, c(2, 3, 1, 4, 5)]

dados <- list(toto1)
resultado <- vector("list", length(toto1))
fator <- 1
dist <- 100 #25, 50, 75, 100, median(radius.new$radius) km [28.097]

require(parallel)
require(doParallel)
for(i in 1:length(dados)) {
  
  #i = 1
  
  ##Getting the list of species data
  list_data <- ConR:::coord.check(XY = dados[[i]],
                                  listing = TRUE)
  
  cl <- makeCluster(detectCores()-2)
  clusterEvalQ(cl, library(sf))
  clusterEvalQ(cl, library(raster))
  clusterEvalQ(cl, library(ConR))
  clusterEvalQ(cl, library(stars))
  
  clusterExport(cl, list("list_data", "fator", "dist"),envir=environment())
  output <- parLapply(cl, 1:length(list_data), fun= function(x) {
    #output <- parLapply(cl, 1:5, fun= function(x) {
    
    # output <- vector("list", length(list_data)) 
    # for(x in 1:length(list_data)) {
    
    res <- try(ConR:::get.patches(
      list_data[[x]],
      cell_size = 2,
      nbe_rep = 0,
      AOO = as.double(unique(list_data[[x]]$AOO)),
      Resol_sub_pop = as.double(unique(list_data[[x]]$radius)),
      subpop_poly = NULL,
      dist_isolated = as.double(unique(list_data[[x]]$radius)) * fator, #for radius.km
      #dist_isolated = dist,
      proj_type = "cea",
      export_shp = FALSE
    ),TRUE)
    
    if(class(res) == "try-error") {
      res <- NULL
      res <- "ERRO"
      names(res) <- unique(list_data[[x]]$tax)
    }
    return(res)
  })
  
  #   cat(x,"\n")
  #   output[[x]] <-res 
  # }    
  
  stopImplicitCluster()
  stopCluster(cl)
  
  output1 <- do.call("c", output)
  df1 <- data.frame(tax = names(output1), frag = as.character(output1), stringsAsFactors = FALSE) 
  resultado[[i]] <- df1
  #saveRDS(df1, paste0("data/tmp",i,".rds"))
} 

resultados <- do.call(rbind.data.frame, resultado)
resultados <- resultados[order(resultados$tax), ]
saveRDS(resultados, paste0("data/04_severe_frag_", dist, "_km_herb.rds"))
#saveRDS(resultados, paste0("data/04_severe_frag_rad_km_herb.rds"))

##Assessing the impact of dist_isolated
AOO <- readRDS("data/04_AOO_herb.rds")
AOO <- tibble::rownames_to_column(AOO, "tax")
dist25 <- readRDS("data/04_severe_frag_25_km_herb.rds")
dist25$frag <- as.double(dist25$frag)
dist50 <- readRDS("data/04_severe_frag_50_km_herb.rds")
dist50$frag <- as.double(dist50$frag)
dist75 <- readRDS("data/04_severe_frag_75_km_herb.rds")
dist75$frag <- as.double(dist75$frag)
dist100 <- readRDS("data/04_severe_frag_100_km_herb.rds")
dist100$frag <- as.double(dist100$frag)
distrad <- readRDS("data/04_severe_frag_rad_km_herb.rds")
distrad$frag <- as.double(distrad$frag)
subpop <- readRDS("data/04_n_subpop_herb.rds")
subpop <- tibble::rownames_to_column(subpop, "species")
names(subpop) <- c("tax", "SUB1")

dist0 <- merge(AOO[, c("tax", "AOO")], subpop[, c("tax", "SUB1")],
               by = "tax", all = TRUE, sort = FALSE) 
table(dist0$tax == dist25$tax)

dist0$frag <- round(100 * dist0$SUB1/(dist0$AOO/4), 2)

jpeg(filename = "figures/Figure_XX.jpg", width = 2250, height = 2000, units = "px", pointsize = 12,
     res = 300, family = "sans", type="cairo", bg="white", quality = 100)

par(mar=c(3,3,1,1), mgp=c(1.5,0.25,0),tcl=-0.2,las=1)
plot(0:6, c(0.8,16,33, 50, 66, 83,100), 
     col= "white", xaxt = "n", cex.lab = 1.2,
     xlab = "Largest distance between subpopulations", ylab="Fragmentation level (%)")#, log = "y")
axis(1, at=c(0.5,1.5,2.5,3.5,4.5,5.5), labels = c("0 km", "25 km","50 km", "75 km", "100 km", "radius"))

#median is represented by the notch
boxplot(as.double(dist0$frag), add=TRUE, notch=TRUE, at=0.5, yaxt = "n")
boxplot(as.double(dist25$frag), add=TRUE, notch=TRUE, at=1.5, yaxt = "n")
boxplot(as.double(dist50$frag), add=TRUE, notch=TRUE, at=2.5, yaxt = "n")
boxplot(as.double(dist75$frag), add=TRUE, notch=TRUE, at=3.5, yaxt = "n")
boxplot(as.double(dist100$frag), add=TRUE, notch=TRUE, at=4.5, yaxt = "n")
boxplot(as.double(distrad$frag), add=TRUE, notch=TRUE, at=5.5, yaxt = "n")
abline(h = 50, lty =3)
legenda <- as.character(c(round((100*table(as.double(dist0$frag)>50)/dim(dist0)[1])[2],1),
                          round((100*table(as.double(dist25$frag)>50)/dim(dist25)[1])[2],1),
                          round((100*table(as.double(dist50$frag)>50)/dim(dist50)[1])[2],1),
                          round((100*table(as.double(dist75$frag)>50)/dim(dist75)[1])[2],1),
                          round((100*table(as.double(dist100$frag)>50)/dim(dist100)[1])[2],1),
                          round((100*table(as.double(distrad$frag)>50)/dim(distrad)[1])[2],1)))
xpos <- c(0.5,1.5,2.5,3.5,4.5,5.5)
for(i in 1:length(xpos)) text(xpos[i], 55, legenda[i])
dev.off()


#### MERGING ALL RESULTS AND SAVING CRITERIA B PARAMETERS ####
rm(list=ls())


#Reading the saved metrics
EOO <- readRDS("data/04_EOO_numeric_herb.rds")
EOO <- tibble::rownames_to_column(EOO, "tax")

AOO <- readRDS("data/04_AOO_herb.rds")
AOO <- tibble::rownames_to_column(AOO, "tax")

subpops <- readRDS("data/04_n_subpop_herb.rds")
subpops <- tibble::rownames_to_column(subpops, "tax")

localities <- readRDS("data/04_number_of_locations_herb.rds")
colnames(localities) <- c("tax", "InsideUC", "OutsideUC")

sev.frag <- merge(AOO, subpops,
                  by = "tax", all = TRUE, sort = FALSE) 
sev.frag$frag<- round(100 * sev.frag$Number_subpop/(sev.frag$AOO/4), 2)

#Checking the species order
table(EOO$tax == AOO$tax)
table(EOO$tax == subpops$tax)
table(EOO$tax == localities$tax)
table(EOO$tax == sev.frag$tax)

####Does cat_criterion_b() understands NA?### TESTING
#EOO = c(200, 100, 400, NA, NA, NA, NA, 100)
#AOO = c(15, 20, 30, 50, 60, 10, 60, 70)
#locations = c(5, 6, 7, 8, 9, 10, 11, 12)
#out <- cat_criterion_b(EOO = EOO, AOO = AOO, locations = locations)
#results_out <- do.call(cbind.data.frame, c(out, stringsAsFactors = FALSE))

#Yes, it does! B1a returns NA also.

#Building the data frame for each level of taxonomic confidence
critB_opt <- cbind.data.frame(EOO[, c("tax", "EOO")],
                              AOO = as.double(AOO[, "AOO"]),
                              Nbe_subPop = as.double(subpops[, "Number_subpop"]),
                              Nbe_loc = as.double(localities[, "OutsideUC"]),
                              Nbe_loc_PA = as.double(localities[, "InsideUC"]),
                              Frag.level = as.double(sev.frag[, "frag"]),
                              stringsAsFactors = FALSE)

saveRDS(critB_opt, "data/04_critB_herb.rds")


##### APPLYING CRITERIA B ##### Need data from criterion C which needs data from pop.decline

#Getting the estimates for criterion B
critB <- readRDS("data/04_critB_herb.rds")
critB$tax <- as.character(critB$tax)

#Getting estimates of species continuing decline based on the index of abundance (Criterion C)
#Usarei o continuing decline do Renato pois ele tem dados para isso que não temos aqui, por conta
#da área de estudo maior — melhor um proxy da AF do que nada
critC <- readRDS("data/criterionC_all_prop_mature.rds")

est.decline <- sapply(strsplit(critC$cont.decline,"\\|"), tail, 1)  
est.decline <- gsub("\\(|\\)|[0-9]", "", est.decline)
est.decline <- gsub(" -", "", est.decline)
table(est.decline, critC$any.decline, useNA = "always")
critC$declineC <- critC$any.decline
critC$declineC[!critC$declineC %in% "Decreasing" & est.decline %in% "Decreasing"] <- "Decreasing" 
critC$declineC[critC$declineC %in% "Increasing" | critC$declineC %in% "Stable"] <- "Not Decreasing" 

#Getting estimates of species continuing decline based on habitat loss (Criterion B)
eoo.decline <- readRDS("rds_files/EOO_hab_loss_2000_2015.rds")
eoo.decline$declineB <- eoo.decline$rel.loss
eoo.decline$declineB[!is.na(eoo.decline$declineB) & eoo.decline$declineB >= 1] <- 1
eoo.decline$declineB[!is.na(eoo.decline$declineB) & as.double(eoo.decline$declineB) < 1 & as.double(eoo.decline$declineB) >= 0] <- 0
eoo.decline$declineB[is.nan(eoo.decline$declineB)] <- NA
eoo.decline$declineB[eoo.decline$declineB %in% 1] <- "Decreasing"
eoo.decline$declineB[eoo.decline$declineB %in% 0] <- "Not Decreasing"
table(eoo.decline$declineB, useNA = "always")

#Correcting the decrease for pioneer species
hab <- read.csv("treeco/threat_habitats.csv", as.is = TRUE)
tmp <- merge(eoo.decline, hab, by.x= "tax", by.y = "Name_submitted", all.x = TRUE, sort = FALSE)
tmp <- tmp[order(tmp$tax),]
table(tmp$tax == eoo.decline$tax)
eoo.decline$net.loss <- eoo.decline$recover - eoo.decline$loss 
hist(eoo.decline$net.loss, nclass = 80)
eoo.decline$declineB[eoo.decline$net.loss >=0.1 & tmp$ecol.group %in% "pioneer"] <- "Not Decreasing"
table(eoo.decline$declineB, useNA = "always")

#Merging any decline info with the criterion B metrics
critB <- merge(critB, critC[,c("species","declineC")],
                   by.x = "tax", by.y = "species", all.x = TRUE, sort = FALSE)
critB <- merge(critB, eoo.decline[,c("tax","declineB")],
                   by = "tax", all.x = TRUE, sort = FALSE)
critB <- critB[order(critB$tax), ]

#Comparing the decline from criterion C and B and correcting if necessary (priority for abundance decline over EOO decline) 
#For assuming, that all other species (rarer) are decreasing as well
hab <- read.csv("treeco/threat_habitats.csv", as.is = TRUE)
tmp <- merge(critB, hab, by.x= "tax", by.y = "Name_submitted", all.x = TRUE, sort = FALSE)
tmp <- tmp[order(tmp$tax), ]
table(critB$tax == tmp$tax)
critB$declineB[critB$declineB %in% "Not Decreasing" & critB$declineC %in%
                     "Decreasing" & !tmp$ecol.group %in% "pioneer"] <- "Decreasing"
table(critB$declineB, critB$declineC, useNA = "always")

#Since some NAs exist, assuming that NA = not decreasing (parsimonious)
critB[c("declineB", "declineC")][is.na(critB[c("declineB", "declineC")])] <- "Not Decreasing"

#Combining the info on number of localities and % in PAs
critB <- 
  critB %>% 
  as_tibble() %>% 
  mutate(nbe_loc_total = Nbe_loc + Nbe_loc_PA) %>% 
  mutate(protected = Nbe_loc_PA/nbe_loc_total*100)

#Creating the severely fragmented column
critB$sever.frag <- 100 * critB$Nbe_subPop/ (critB$AOO / 4) > 50


#### PERFOMING THE ASSESSMENTES OF CRITERION B ####
results_Cb <- cat_criterion_b(EOO = critB$EOO,
                              AOO = critB$AOO,
                              locations = critB$nbe_loc_total,
                              sever.frag = critB$sever.frag,
                              protected = critB$Nbe_loc_PA, 
                              decline = critB$declineB,
                              protected.threshold = 100 
)
                                 
sum((100 * table(results_Cb$ranks_B, useNA = "always")/dim(critB)[1])[c(1,2,4)]) #31.67% #26.84% (herb)

#Saving the results
results_Cb <- do.call(cbind.data.frame, c(results_Cb, stringsAsFactors = FALSE))
critB_opt.all <- critB[,c("tax","EOO","AOO","Nbe_subPop","nbe_loc_total","protected","declineB","sever.frag")]
critB_opt.all[, c("category_B", "category_B_code","category_B1","category_B2")] <- NA_character_
critB_opt.all[, c("category_B", "category_B_code","category_B1","category_B2")] <- results_Cb
critB_opt.all[is.na(critB_opt.all$AOO),]

saveRDS(critB_opt.all, "data/critB_herb.rds")
saveRDS(as.data.frame(critB_opt.all), "data/criterionB_herb.rds")

rm(list=ls())
