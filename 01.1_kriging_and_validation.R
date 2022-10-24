rm(list = ls())

{
library(tidyverse)
library(performance)
library(gstat)
library(sp)
library(automap)
library(raster)
library(rgdal)
library(sf)
library(MASS)
library(lattice)
}

#Loading data 
grid_env = readRDS('data_old_stuff/00_grid_env.rds')
points = readRDS('data/00_points_filtered_nm.rds')
sc_buffer = readRDS('data_old_stuff/00_sc_buffer.rds')
sc_wo_grassland = readRDS('data_old_stuff/00_sc_wo_grassland.rds')

#Selecting only chosen predictors from best subset selection script
points = points[, c("SiteCode.x", "DA", "DA_new", "temp", "Tmin.wclim", 
                    "um.rel", "CWD", "EnvStress", "long1", "lat1")]

#Removing outliers DA/ha
points = points %>% filter(DA_new <= 1200)


##### INTERPOLATING TREE DENSITY/HECTARE #####

#Variography using automap library
coordinates(points) = ~long1+lat1
crs(points) = crs(sc_buffer)

#Variables already chosen using best subset selection (see script)
v <- autofitVariogram(DA_new ~ temp + Tmin.wclim, points)
jpeg(filename = "figures/variogram.jpg", width = 2500, height = 2250, units = "px", pointsize = 12,
     res = 300, family = "sans", type = "cairo", bg = "white", quality = 100)
plot(v)
dev.off()
v$var_model

#Universal Kriging
crs(grid_env) = crs(points)

gstat_krig <- gstat(formula = DA_new ~ temp + Tmin.wclim, 
                    locations = points, model = v$var_model)

grid_env = as.data.frame(grid_env)
grid_env = grid_env[, c("temp", "Min.Temp.Coldest", "Longitude", "Latitude")]
coordinates(grid_env) = ~Longitude + Latitude
gridded(grid_env) = TRUE
crs(grid_env) = crs(points)
colnames(grid_env@data) = c("temp", "Tmin.wclim")
k_output <- predict(gstat_krig, grid_env)
k_clipped <- crop(k_output, sc_wo_grassland)
plot(k_clipped, main = "DA = 1200")

#Adding upper and lower confidence interval from the mean
k_clipped$lower = (k_clipped$var1.pred - 1.96*(sqrt(k_clipped$var1.pred)/sqrt(dim(points)[1])))
k_clipped$upper = (k_clipped$var1.pred + 1.96*(sqrt(k_clipped$var1.pred)/sqrt(dim(points)[1])))
k_raster <- raster(k_clipped)
saveRDS(k_clipped, 'data/01_kriging_output.rds')
writeRaster(k_raster, 'Resultados/kriging_output.tif', overwrite = TRUE)
writeOGR(k_clipped, 'Resultados/kriging_output.shp', driver = 'ESRI Shapefile', layer = 'var1.pred')

#Adding upper and lower values to forest cover data
df <- as.data.frame(sc_wo_grassland)
df <- df[, -c(7, 8)]
df <- df[, 1:41]
df$krig_lower <- k_clipped$lower
df$krig_fit <- k_clipped$var1.pred
df$krig_upper <- k_clipped$upper
df <- df[, c(1:6, 42:44, 7:41)]
write.csv(df, 'sc_wo_grassland_df.csv')

#Testing how many trees we had in SC (1850)
cellStats(k_raster, sum)
1873895 * 2500 #original forest cover (25km² all covered by forest: 4.68 billions)

##### Kriging with cross validation (LOOCV) to validate the DA/ha model #####
k.cv <- krige.cv(DA_new ~ temp + Tmin.wclim, points, grid_env, model=v$var_model)

# mean error, ideally 0:
mean(k.cv$residual)
# MSPE, ideally small
mean(k.cv$residual^2)
# Mean square normalized error, ideally close to 1
mean(k.cv$zscore^2)
# correlation observed and predicted, ideally 1
cor(k.cv$observed, k.cv$observed - k.cv$residual)
# correlation predicted and residual, ideally 0
cor(k.cv$observed - k.cv$residual, k.cv$residual)
# }

coordinates(k.cv) = ~long1+lat1

#Linear model between observed vs predicted
model_lm <- lm(observed ~ var1.pred, data = k.cv)

#Validating the model using RMSE and adj.R-squared #Good enough model for DA/ha
model_performance(model_lm)
check_model(model_lm)

#Plotting the linear model
xv1 <- data.frame(k.cv$observed, k.cv$var1.pred)
xv1 <- round(xv1, 0)

#add R-squared and/or RMSE
jpeg(filename = "Resultados/linear_model.jpeg", width = 2500, height = 2150, units = "px", pointsize = 12,
res = 300, family = "sans", type="cairo", bg="white", quality = 100)
ggplot(data = xv1, aes(x = k.cv.observed, y = k.cv.var1.pred)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
labs(x = "Observed DA/ha", y = "Predicted DA/ha", title = "Linear regression model for predicted vs observed DA/ha") +
  theme_light(base_size = 16) +
  theme(plot.title = element_text(hjust=0.5)) +
  annotate("text", x = 1000, y = 450, label = "italic(R)^2 == 0.16", parse = TRUE, size = 6)
dev.off()

##### End #####
rm(list = ls())
