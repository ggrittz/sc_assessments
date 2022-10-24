rm(list=ls())

library(dplyr)
library(leaps)

data = readRDS('data/00_points_filtered_nm.rds')
data = data[, c("DA_new", 'temp', 'Temp.Seasonality', 'Tmin.wclim', 'ppt', 'um.rel', 
                'EnvStress', 'DirRad_Eq2.med', 'DECLIV', 'rad.wclim', 'HumanInfluence')]

#Are all values numeric?
sapply(data, is.factor)
#data = data %>% mutate_at(vars(CWD, EnvStress, DECLIV), as.character)
#data = data %>% mutate_at(vars(CWD, EnvStress, DECLIV), as.numeric)
#sum(is.na(data))

data[,2:11] <- scale(data[, 2:11])

#Generating the models
regfit.full = regsubsets(DA_new ~ temp+Temp.Seasonality+Tmin.wclim+ppt+um.rel+
                         EnvStress+DirRad_Eq2.med+DECLIV+rad.wclim+HumanInfluence, data)
summary(regfit.full)

#preparing summary model selections methods
reg.summary = summary(regfit.full)
names(reg.summary)

#Model selection using Cp and BIC: which.min() function indicates the smallest value for each metric
#Mallow's CP
plot(reg.summary$cp, xlab = 'Number of Variables', ylab = "Mallow's Cp", type = 'l')
which.min(reg.summary$cp)
points(3, reg.summary$cp[3], col = 'red', cex = 2, pch = 20)

#BIC
which.min(reg.summary$bic)
svg('Resultados/subset_selection.svg', width = 10, height = 8)
plot(reg.summary$bic, xlab = 'Number of Variables', ylab = 'BIC', type = 'l', main = "Best subset selection of predictors")
points(2, reg.summary$bic[2], col = 'red', cex = 2, pch = 20)
dev.off()

#built-in function of regsubsets() to display selected variables for the best model
plot(regfit.full, scale='bic', ylab = "BIC")

#Which variables were selected? #temp, Tmin.wclim
coef(regfit.full, 2)
