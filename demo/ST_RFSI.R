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
data(dtempc_ogimet)

dtempc <- dtempc[complete.cases(dtempc),]
summary(dtempc)
dtempc.sf <- st_as_sf(dtempc, coords = c("lon", "lat"), crs = 4326, agr = "constant")
# dtempc.sf <- st_transform(dtempc.sf, 32634)
dtempc.sftime <- st_sftime(dtempc.sf)
class(dtempc.sftime)
st_crs(dtempc.sftime)
head(unique(st_time(dtempc.sftime)))
# plot(dtempc.sftime[, "tmean"])

#################### rfsi ####################

# fm.RFSI <- as.formula("om ~ dist + soil") # covariates without spatial covariates (obs1, dist1, ...)
fm.RFSI <- as.formula("tmean ~ gtt + dem + twi + cdate + doy")
data <- dtempc.sftime
# data <- cbind(st_drop_geometry(dtempc.sftime), st_coordinates(dtempc.sftime[, "geometry"]))
# data = dtempc.sf
# data.staid.x.y.z <- c("staid","X","Y","time")

rfsi_model <- rfsi(formula = fm.RFSI,
                   data = data, # dtempc
                   # data.staid.x.y.z = data.staid.x.y.z, # only if class(data) == data.frame
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
data.staid.x.y.z <- c("staid","X","Y","time")
class(data)

# Tune an RFSI model
rfsi_tuned <- tune.rfsi(formula = fm.RFSI,
                        data = data, # meuse.df (use data.staid.x.y.z)
                        data.staid.x.y.z = data.staid.x.y.z, # only if class(data) == data.frame
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

# prepare covariates
data(dem_twi_srb)
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

# Cross-validation of RFSI
rfsi_cv <- cv.rfsi(formula=fm.RFSI, # without nearest obs
                   data = data, # meuse.df (use data.staid.x.y.z)
                   # data.staid.x.y.z = c("id", "x", "y", NA), # only if class(data) == data.frame
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




