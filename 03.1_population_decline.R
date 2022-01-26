#### POPULATION DECLINE FOR MISSING YEARS ####
library(ConR)
library(tidyverse)

#devtools::install_github("gdauby/ConR@devel")

#### PREPARING DATA ####
#lower conf interval
lower_popsize <- readRDS("data/03_timeseries_lower_popsize.rds")
lower_popsize <- tibble::rownames_to_column(lower_popsize, "species")

#fit (mean) conf interval
fit_popsize <- readRDS("data/03_timeseries_fit_popsize.rds")
fit_popsize <- tibble::rownames_to_column(fit_popsize, "species")
fit_popsize$species <- gsub("\\.", " ", fit_popsize$species)

#upper conf interval
upper_popsize <- readRDS("data/03_timeseries_upper_popsize.rds")
upper_popsize <- tibble::rownames_to_column(upper_popsize, "species")

####POPULATION SIZE####
colClean <- function(x){ colnames(x) <- gsub("trees_", "", colnames(x)); x} 
lower_popsize = colClean(lower_popsize)
fit_popsize = colClean(fit_popsize)
upper_popsize = colClean(upper_popsize)

known.years = c(1850, 1985:2019)

## OBTAINING POPULATION SIZES FOR MISSING YEARS ##
## Getting pop. sizes for all possible generation lengths
gen.lengths <- c(10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100)
pos.gen.length <- sort(unique(gen.lengths * rep(c(1:3), each=length(gen.lengths))), decreasing = TRUE)
miss.years <- 2019 - pos.gen.length
miss.years1 <- miss.years[miss.years < 1850]
miss.years2 <- miss.years[miss.years > 1850 & miss.years < 1985]
miss.years3 <- c(miss.years1, miss.years2)

#### Looping the pop. decline trends for each species and time interval ####
require(snow)
require(doSNOW)
require(foreach)

dados <- list(lower_popsize[, 2:37], fit_popsize[, 2:37], upper_popsize[, 2:37])
results <- vector("list", length(dados))

for(i in 1:length(dados)) {
  #defining the object for the loop
  pop_data <- dados[[i]]
  
  #Setting the loop parameters
  cl <- snow::makeSOCKcluster(4)
  doSNOW::registerDoSNOW(cl)
  `%d%` <- foreach::`%dopar%`
  x <- NULL
  output <-
    foreach::foreach(
      x = 1:dim(pop_data)[1],
      #.combine = 'c',
      .options.snow = NULL
    ) %d% {
      
      res1 <- ConR::pop.decline.fit(pop.size = pop_data[x, ],
                                     years = colnames(dados[[1]]),
                                     models = "all", 
                                     project.years = miss.years3,
                                     plot.fit = FALSE)
      
      res1$data$Modelo <- attributes(res1$best.model)$best.model.name
      dados <- rbind.data.frame(res1$data, stringsAsFactors = FALSE) 
    }
  snow::stopCluster(cl)
  
  #Saving the predictions for each table (mean, lower and upper) 
  results[[i]] <- output 
}

#Setting species names
for(i in 1:length(results)) names(results[[i]]) <- fit_popsize$species

## Putting data in the ConR format
years <- results[[1]][[1]]$Year
ncols <- length(years)
spp <- names(results[[1]])
nrows <- length(spp)
high.pop.conR <- mean.pop.conR <- low.pop.conR <- matrix(NA, ncol = ncols, nrow = nrows,
                                                         dimnames = list(spp, years))
for(x in 1:length(results[[1]])) {
  low.pop.conR[x, ] <- results[[1]][[x]]$Predicted
  mean.pop.conR[x, ] <- results[[2]][[x]]$Predicted
  high.pop.conR[x, ] <- results[[3]][[x]]$Predicted
}

#Converting into data.frames
low.pop.conR <- cbind.data.frame(species = rownames(low.pop.conR), low.pop.conR, 
                                 row.names = NULL, stringsAsFactors = FALSE)

mean.pop.conR <- cbind.data.frame(species = rownames(mean.pop.conR), mean.pop.conR, 
                                  row.names = NULL, stringsAsFactors = FALSE)

high.pop.conR <- cbind.data.frame(species = rownames(high.pop.conR), high.pop.conR, 
                                  row.names = NULL, stringsAsFactors = FALSE)

#### SAVING ###

#Saving the estimated and infered populations (BOTH "OBSERVED" AND ESTIMATED/INTERPOLATED POP SIZES)
saveRDS(results[[1]], "data/threat_low_pop_sizes_infer.rds")
saveRDS(results[[2]], "data/threat_mean_pop_sizes_infer.rds")
saveRDS(results[[3]], "data/threat_high_pop_sizes_infer.rds")

#Saving the estimated and infered populations in the ConR format
saveRDS(low.pop.conR, "data/threat_low_pop_sizes_for_ConR.rds")
saveRDS(mean.pop.conR, "data/threat_mean_pop_sizes_for_ConR.rds")
saveRDS(high.pop.conR, "data/threat_high_pop_sizes_for_ConR.rds")

rm(list=ls())
