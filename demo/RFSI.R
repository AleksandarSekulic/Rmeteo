library(meteo)
library(sp)
library(sf)
library(sftime)
library(terra)
library(gstat)
library(plyr)
library(xts)
library(snowfall)
library(doParallel)
library(CAST)
library(ranger)

# preparing data
demo(meuse, echo=FALSE)
meuse <- meuse[complete.cases(meuse@data),]

#################### near.obs ####################
locations = st_as_sf(meuse, coords = c("x", "y"), crs = 28992, agr = "constant") # sf

nearest_obs <- near.obs(locations = locations, # from
                        locations.x.y = c("x","y"),
                        observations = locations, # to
                        observations.x.y=c("x","y"),
                        obs.col = "zinc",
                        n.obs = 10, # number of nearest observations
                        rm.dupl = TRUE) # if TRUE, the observations themselves will be excluded from spatial covariates 
str(nearest_obs)
summary(nearest_obs)

#################### rfsi ####################

# fm.RFSI <- as.formula("om ~ dist + soil") # covariates without spatial covariates (obs1, dist1, ...)
fm.RFSI <- as.formula("zinc ~ dist + soil + ffreq")
data = st_as_sf(meuse, coords = c("x", "y"), crs = 28992, agr = "constant")
# data = terra::vect(meuse)
# data <- as.data.frame(meuse)
# data$id = 1:nrow(data)
# data.staid.x.y.z <- c("id","x","y",NA)
class(data)

rfsi_model <- rfsi(formula = fm.RFSI,
                   data = data,
                   # data.staid.x.y.z = data.staid.x.y.z, # only if class(data) == data.frame
                   n.obs = 5, # number of nearest observations
                   # s.crs = st_crs(data), # nedded only if the coordinates are lon/lat (WGS84)
                   # p.crs = st_crs(data), # nedded only if the coordinates are lon/lat (WGS84)
                   cpus = detectCores()-1,
                   progress = TRUE,
                   # ranger parameters
                   importance = "impurity",
                   seed = 42,
                   num.trees = 250,
                   mtry = 5,
                   splitrule = "variance",
                   min.node.size = 5,
                   sample.fraction = 0.95,
                   quantreg = FALSE) # quantile regression model

rfsi_model
# OOB prediction error (MSE):       47758.14 
# R squared (OOB):                  0.6435869 
# Note that OOB error statistics are biased and should not be considered as accuracy metrics (they do not show spatial accuracy)!
# The proper way to assess accuaracy of the RFSI model is by using the nested k-fold cross-validation (cv.rfsi function)
sort(rfsi_model$variable.importance)

#################### pred.rfsi ####################

# newdata <- as.data.frame(meuse.grid)
# newdata$id <- 1:nrow(newdata)
# newdata <- terra::vect(meuse.grid)
newdata <- terra::rast(meuse.grid)
class(newdata)

rfsi_prediction <- pred.rfsi(model = rfsi_model,
                             data = data,
                             obs.col = "zinc",
                             # data.staid.x.y.z = data.staid.x.y.z, # only if class(data) == data.frame
                             newdata = newdata, # meuse.grid.df (use newdata.staid.x.y.z)
                             # newdata.staid.x.y.z = c("id", "x", "y", NA), # only if class(newdata) == data.frame
                             output.format = "SpatRaster", # "sf", # "SpatVector", 
                             # s.crs = st_crs(data), # NA
                             # newdata.s.crs = st_crs(data), # NA
                             # p.crs = st_crs(data), # NA
                             cpus = 1, # detectCores()-1,
                             progress = TRUE,
                             soil3d=FALSE
                             # type = "quantiles",
                             # quantiles = c(0.1, 0.5, 0.9)
)

class(rfsi_prediction)
names(rfsi_prediction)
summary(rfsi_prediction)

plot(rfsi_prediction)
# plot(rfsi_prediction['pred'])
# plot(rfsi_prediction['quantile..0.1'])
# plot(rfsi_prediction['quantile..0.5'])
# plot(rfsi_prediction['quantile..0.9'])

#################### tune.rfsi ####################

# making tgrid
n.obs <- 1:6
min.node.size <- 2:10
sample.fraction <- seq(1, 0.632, -0.05) # 0.632 without / 1 with replacement
splitrule <- "variance"
ntree <- 250 # 500
mtry <- 3:(2+2*max(n.obs))
tgrid = expand.grid(min.node.size=min.node.size, num.trees=ntree,
                    mtry=mtry, n.obs=n.obs, sample.fraction=sample.fraction)

# Tune an RFSI model
rfsi_tuned <- tune.rfsi(formula = fm.RFSI,
                        data = data,
                        # data.staid.x.y.z = data.staid.x.y.z, # only if class(data) == data.frame
                        # s.crs = NA,
                        # p.crs = NA,
                        tgrid = tgrid, # combinations for tuning
                        tgrid.n = 20, # number of randomly selected combinations from tgrid for tuning
                        tune.type = "LLO", # Leave-Location-Out CV
                        k = 5, # number of folds
                        acc.metric = "RMSE", # R2, CCC, MAE
                        fit.final.model = TRUE,
                        cpus = detectCores()-1,
                        progress = TRUE,
                        # ranger parameters
                        importance = "impurity",
                        seed = 42
                        # type = "quantiles",
                        # quantiles = c(0.1, 0.5, 0.9)
)

rfsi_tuned$combinations
rfsi_tuned$tuned.parameters
# min.node.size num.trees mtry n.obs sample.fraction     RMSE
# 3701             3       250    6     5            0.75 222.6752
rfsi_tuned$final.model
# OOB prediction error (MSE):       46666.51 
# R squared (OOB):                  0.6517336 

#################### cv.rfsi ####################

# Cross-validation of RFSI
rfsi_cv <- cv.rfsi(formula=fm.RFSI, # without nearest obs
                   data = data,
                   # data.staid.x.y.z = c("id", "x", "y", NA), # only if class(data) == data.frame
                   # s.crs=NA,
                   # p.crs=NA,
                   tgrid = tgrid, # combinations for tuning
                   tgrid.n = 5, # number of randomly selected combinations from tgrid for tuning
                   tune.type = "LLO", # Leave-Location-Out CV
                   k = 5, # number of outer and inner folds
                   seed = 42,
                   acc.metric = "RMSE", # R2, CCC, MAE
                   output.format = "sf", # "SpatVector", # "data.frame",
                   cpus=detectCores()-1,
                   progress=1,
                   # ranger parameters
                   importance = "impurity")

summary(rfsi_cv)
rfsi_cv$dif <- rfsi_cv$obs - rfsi_cv$pred
plot(rfsi_cv["dif"])
acc.metric.fun(rfsi_cv$obs, rfsi_cv$pred, "R2") # 0.5948994
acc.metric.fun(rfsi_cv$obs, rfsi_cv$pred, "RMSE") # 232.2175
acc.metric.fun(rfsi_cv$obs, rfsi_cv$pred, "MAE") # 0.3328514
acc.metric.fun(rfsi_cv$obs, rfsi_cv$pred, "CCC") # 0.7484076

