library(meteo)
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
# load(file = "/media/sekulic/e42d0594-47d8-4b72-99f7-bf3809dbeb4f/__DOKTORAT/MeteoSerbia1km/ogimet/ogimet_serbia08_tmean.rda")
# dtempc <- ogimet_serbia[, c(1:12)]
# dtempc <- dtempc[as.numeric(format(dtempc$date, "%Y")) == 2019, ]
# names(dtempc)[6] <- "time"
# saveRDS(dtempc, file="/home/sekulic/Projects/Rmeteo/data/dtempc.rds")
# dtempc <- readRDS("/home/sekulic/Projects/Rmeteo/data/dtempc.rds")
# save(dtempc, file = "/home/sekulic/Projects/Rmeteo/data/dtempc.rda")
# stations <- unique(dtempc[, c(1,4,3,5,2,8,9)])
# save(stations, file = "/home/sekulic/Projects/Rmeteo/data/stations.rda")
data(dtempc)
# data(stations)

dtempc <- dtempc[complete.cases(dtempc),]
# dtempc <- join(dtempc, stations)
summary(dtempc)
dtempc.sf <- st_as_sf(dtempc, coords = c("lon", "lat"), crs = 4326, agr = "constant")
# dtempc.sf <- st_transform(dtempc.sf, 32634)
dtempc.sftime <- st_sftime(dtempc.sf)
class(dtempc.sftime)
st_crs(dtempc.sftime)
head(unique(st_time(dtempc.sftime)))
# plot(dtempc.sftime[, "tmean"])

# serbia.bbox <- st_bbox(c(xmin = 18, xmax = 22.5, ymax = 40, ymin = 46), crs = st_crs(4326))
# serbia.sftime <- st_filter(dtempc.sftime, st_as_sfc(serbia.bbox))
# class(serbia.sftime)
# plot(serbia.sftime)

#################### rfsi ####################
# source("~/Projects/Rmeteo/R/data.prepare.R")
# source("~/Projects/Rmeteo/R/near.obs.R")
# source("~/Projects/Rmeteo/R/rfsi.R")

# fm.RFSI <- as.formula("om ~ dist + soil") # covariates without spatial covariates (obs1, dist1, ...)
fm.RFSI <- as.formula("tmean ~ gtt + dem + twi + cdate + doy")
data <- dtempc.sftime
# data <- cbind(st_drop_geometry(dtempc.sftime), st_coordinates(dtempc.sftime[, "geometry"]))
# data = dtempc.sf
# data.staid.x.y.z <- c("staid","X","Y","time")

rfsi_model <- rfsi(formula = fm.RFSI,
                   data = data, # dtempc
                   # data.staid.x.y.z = data.staid.x.y.z, # only if class(data) == data.frame
                   zero.tol = 0,
                   n.obs = 5, # number of nearest observations
                   s.crs = st_crs(4326), # nedded only if the coordinates are lon/lat (WGS84)
                   p.crs = st_crs(32634), # nedded only if the coordinates are lon/lat (WGS84)
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
# OOB prediction error (MSE):       0.8872976 
# R squared (OOB):                  0.9872141 
# Note that OOB error statistics are biased and should not be considered as accuracy metrics (they do not show spatial accuracy)!
# The proper way to assess accuaracy of the RFSI model is by using the nested k-fold cross-validation (cv.rfsi function)
sort(rfsi_model$variable.importance)

#################### tune.rfsi ####################
# source("~/Projects/Rmeteo/R/near.obs.R")
# source("~/Projects/Rmeteo/R/data.prepare.R")
# source("~/Projects/Rmeteo/R/rfsi.R")
# source("~/Projects/Rmeteo/R/pred.rfsi.R")
# source("~/Projects/Rmeteo/R/tune.rfsi.R")
# source("~/Projects/Rmeteo/R/acc.metric.fun.R")

# making tgrid
n.obs <- 1:10
min.node.size <- 2:10
sample.fraction <- seq(1, 0.632, -0.05) # 0.632 without / 1 with replacement
splitrule <- "variance"
ntree <- 250 # 500
mtry <- 3:(2+2*max(n.obs))
tgrid = expand.grid(min.node.size=min.node.size, num.trees=ntree,
                    mtry=mtry, n.obs=n.obs, sample.fraction=sample.fraction)

fm.RFSI <- as.formula("tmean ~ gtt + dem + twi + cdate + doy")
data <- dtempc.sftime
# data <- cbind(st_drop_geometry(dtempc.sftime), st_coordinates(dtempc.sftime[, "geometry"]))
# data = dtempc.sf
# data.staid.x.y.z <- c("staid","X","Y","time")
class(data)

# Tune an RFSI model
rfsi_tuned <- tune.rfsi(formula = fm.RFSI,
                        data = data, # meuse.df (use data.staid.x.y.z)
                        data.staid.x.y.z = data.staid.x.y.z, # only if class(data) == data.frame
                        zero.tol = 0,
                        s.crs = st_crs(4326), # NA
                        p.crs = st_crs(32634), # NA
                        tgrid = tgrid, # combinations for tuning
                        tgrid.n = 10, # number of randomly selected combinations from tgrid for tuning
                        tune.type = "LLO", # Leave-Location-Out CV
                        k = 5, # number of folds
                        acc.metric = "RMSE", # R2, CCC, MAE
                        fit.final.model = TRUE,
                        cpus = detectCores()-1,
                        progress = TRUE,
                        # ranger parameters
                        importance = "impurity",
                        seed = 42)
# type = "quantiles",
# quantiles = c(0.1, 0.5, 0.9)

rfsi_tuned$combinations
rfsi_tuned$tuned.parameters
# min.node.size num.trees mtry n.obs sample.fraction     RMSE
# 1165             5       250   12     7               1 1.531365
rfsi_tuned$final.model
# OOB prediction error (MSE):       0.8927979 
# R squared (OOB):                  0.9871348

#################### pred.rfsi ####################
source("~/Projects/Rmeteo/R/tgeom2STFDF.R")
source("~/Projects/Rmeteo/R/pred.rfsi.R")

# prepare covariates
# data(regdata)
data(dem_twi_srb)
# dem_twi_srb <- readRDS("/home/sekulic/Projects/Rmeteo/data/dem_twi_srb.rds")
# save(dem_twi_srb, file="/home/sekulic/Projects/Rmeteo/data/dem_twi_srb.rda")
# data(tregcoef)
# dem <- rast("/media/sekulic/e42d0594-47d8-4b72-99f7-bf3809dbeb4f/__DOKTORAT/MeteoSerbia1km/dem_twi/dem.tif")
# twi <- rast("/media/sekulic/e42d0594-47d8-4b72-99f7-bf3809dbeb4f/__DOKTORAT/MeteoSerbia1km/dem_twi/twi.tif")
# dem.twi <- c(dem, twi)
# plot(dem.twi["dem"])
# saveRDS(dem.twi, file="/home/sekulic/Projects/Rmeteo/data/dem.twi.rds")
# terra::writeRaster(regdata, "/home/sekulic/Projects/Rmeteo/data/regdata.tif", filetype = "GTiff", overwrite = TRUE)
# regdata <- rast("/home/sekulic/Projects/Rmeteo/data/regdata.tif")
# covariates <- rast(readRDS("/home/sekulic/Projects/Rmeteo/data/dem_twi_srb.rds"))
covariates <- rast(dem_twi_srb)
date = "2019-01-01"
doy = as.integer(strftime(as.POSIXct(paste(date), format="%Y-%m-%d"), format = "%j"))
cdate = floor(unclass(as.POSIXct(as.POSIXct(paste(date), format="%Y-%m-%d")))/86400)
covariates$doy = doy
covariates$cdate = cdate
gtt <- temp_geom(day=doy,
                 fi = crds(covariates, na.rm=FALSE)[,2],
                 variable="mean",
                 ab=NULL)
covariates$gtt = gtt
names(covariates)



# model = rfsi_model
# obs.col = "tmean"
# data <- dtempc.sftime
# data.staid.x.y.z = NULL
# z.value = date
# output.format = "SpatRaster"
# zero.tol = 0
# s.crs = st_crs(4326)
# newdata.s.crs = st_crs(4326)
# p.crs = st_crs(32634)
# cpus = 1
# progress = TRUE
# soil3d=FALSE

newdata <- covariates
# newdata.staid.x.y.z = NULL

# newdata <- cbind(crds(covariates, df=T, na.rm=T), as.data.frame(covariates, na.rm=T))
# newdata$id <- 1:nrow(newdata)
# newdata$time <- date
# newdata.staid.x.y.z = c("id", "x", "y", "time")
# newdata <- terra::vect(meuse.grid)
class(newdata)

rfsi_prediction <- pred.rfsi(model = rfsi_tuned$final.model, # rfsi_model,
                             data = data,
                             obs.col = "tmean",
                             # data.staid.x.y.z = data.staid.x.y.z, # only if class(data) == data.frame
                             newdata = newdata, # meuse.grid.df (use newdata.staid.x.y.z)
                             # newdata.staid.x.y.z = newdata.staid.x.y.z, # only if class(newdata) == data.frame
                             z.value = date, # c("2019-01-01", "2019-01-02")
                             output.format = "SpatRaster", # "sf", # "SpatVector", 
                             zero.tol = 0,
                             s.crs = st_crs(4326), # NA
                             newdata.s.crs = st_crs(4326), # NA
                             p.crs = st_crs(32634), # NA
                             cpus = 1, # detectCores()-1,
                             progress = TRUE,
                             soil3d=FALSE
                             # type = "quantiles",
                             # quantiles = c(0.1, 0.5, 0.9)
)

class(rfsi_prediction)
names(rfsi_prediction)
summary(rfsi_prediction)

plot(rfsi_prediction$`2019-01-01`)

#################### cv.rfsi ####################
source("~/Projects/Rmeteo/R/near.obs.R")
source("~/Projects/Rmeteo/R/data.prepare.R")
source("~/Projects/Rmeteo/R/rfsi.R")
source("~/Projects/Rmeteo/R/pred.rfsi.R")
source("~/Projects/Rmeteo/R/tune.rfsi.R")
source("~/Projects/Rmeteo/R/acc.metric.fun.R")
source("~/Projects/Rmeteo/R/cv.rfsi.R")

# Cross-validation of RFSI
rfsi_cv <- cv.rfsi(formula=fm.RFSI, # without nearest obs
                   data = data, # meuse.df (use data.staid.x.y.z)
                   # data.staid.x.y.z = c("id", "x", "y", NA), # only if class(data) == data.frame
                   zero.tol=0,
                   s.crs = st_crs(4326), # NA
                   p.crs = st_crs(32634), # NA
                   tgrid = tgrid, # combinations for tuning
                   tgrid.n = 5, # number of randomly selected combinations from tgrid for tuning
                   tune.type = "LLO", # Leave-Location-Out CV
                   k = 5, # number of folds
                   seed = 42,
                   acc.metric = "RMSE", # R2, CCC, MAE
                   output.format = "SpatVector", # "sf", # "data.frame",
                   cpus=detectCores()-1,
                   progress=1,
                   # ranger parameters
                   importance = "impurity")

summary(rfsi_cv)
rfsi_cv$dif <- rfsi_cv$obs - rfsi_cv$pred
rfsi_cv_st <- aggregate(rfsi_cv, by='staid')
plot(rfsi_cv_st, 'mean_dif')
acc.metric.fun(rfsi_cv$obs, rfsi_cv$pred, "R2") # 0.9652058
acc.metric.fun(rfsi_cv$obs, rfsi_cv$pred, "RMSE") # 1.55386
acc.metric.fun(rfsi_cv$obs, rfsi_cv$pred, "MAE") # 0.08545364
acc.metric.fun(rfsi_cv$obs, rfsi_cv$pred, "CCC") # 0.9822148




