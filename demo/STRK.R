library(meteo)
library(sp)
library(spacetime)
library(gstat)
library(plyr)
library(xts)
library(snowfall)
library(doParallel)
library(CAST)

# preparing data
data(dtempc) 
data(stations)
data(regdata) # covariates, made by mete2STFDF function

regdata@sp@proj4string <- CRS('+proj=longlat +datum=WGS84')
data(tvgms) # ST variogram models
data(tregcoef) # MLR coefficients

lonmin=18 ;lonmax=22.5 ; latmin=40 ;latmax=46
serbia = point.in.polygon(stations$lon, stations$lat, c(lonmin,lonmax,lonmax,lonmin), 
                          c(latmin,latmin,latmax,latmax))
st = stations[ serbia!=0, ] # stations in Serbia approx.
crs = CRS('+proj=longlat +datum=WGS84')

# create STFDF
stfdf <- meteo2STFDF(obs = dtempc,
                     stations = st,
                     crs = crs)

#################### pred.strk ####################

# Calculate prediction of mean temperatures for "2011-07-05" and "2011-07-06" 
# global model is used for regression and variogram

### Example with STFDF and without parallel processing
results <- pred.strk(data = stfdf, # observations
                     newdata = regdata, # prediction locations with covariates
                     # newdata = regdata[,2,drop=FALSE], # for one day only
                     output.format = "STFDF", # data.frame | sf | sftime | SpatVector | SpatRaster 
                     reg.coef = tregcoef[[1]], # MLR coefficients
                     vgm.model = tvgms[[1]], # STRK variogram model
                     sp.nmax = 20,
                     time.nmax = 2,
                     computeVar=TRUE
)
class(results)
# plot prediction
results@sp=as(results@sp,'SpatialPixelsDataFrame')
stplot(results[,,"pred", drop= FALSE], col.regions=bpy.colors())
stplot(results[,,"var", drop= FALSE], col.regions=bpy.colors())

### Example with data.frames and parallel processing - SpatRaster output
# create data.frame
stfdf.df <- join(dtempc, st)
summary(stfdf.df)
regdata.df <- as.data.frame(regdata)

results <- pred.strk(data = stfdf.df,
                     obs.col = 3,
                     data.staid.x.y.z = c(1,4,5,2),
                     newdata = regdata.df,
                     newdata.staid.x.y.z = c(3,1,2,4),
                     crs = CRS("EPSG:4326"),
                     output.format = "SpatRaster", # STFDF |data.frame | sf | sftime | SpatVector
                     reg.coef = tregcoef[[1]],
                     vgm.model = tvgms[[1]],
                     sp.nmax = 20,
                     time.nmax = 2,
                     parallel.processing = TRUE,
                     pp.type = "doParallel", # "snowfall"
                     cpus = detectCores()-1,
                     computeVar = TRUE,
                     progress = TRUE
)

# plot prediction
plot(results$`2011-07-06`[["pred"]])
plot(results$`2011-07-06`[["var"]])

#################### cv.strk ####################

# Cross-validation for mean temperature for days "2011-07-05" and "2011-07-06" 
# global model is used for regression and variogram

# Overlay observations with covariates
time <- index(stfdf@time)
covariates.df <- as.data.frame(regdata)
names_covar <- names(tregcoef[[1]])[-1]
for (covar in names_covar){
  nrowsp <- length(stfdf@sp)
  regdata@sp=as(regdata@sp,'SpatialPixelsDataFrame')
  ov <- sapply(time, function(i) 
    if (covar %in% names(regdata@data)) {
      if (as.Date(i) %in% as.Date(index(regdata@time))) {
        over(stfdf@sp, as(regdata[, i, covar], 'SpatialPixelsDataFrame'))[, covar]
      } else {
        rep(NA, length(stfdf@sp))
      }
    } else {
      over(stfdf@sp, as(regdata@sp[covar], 'SpatialPixelsDataFrame'))[, covar]
    }
  )
  ov <- as.vector(ov)
  if (all(is.na(ov))) {
    stop(paste('There is no overlay of data with covariates!', sep = ""))
  }
  stfdf@data[covar] <- ov
}

# Remove stations out of covariates
for (covar in names_covar){
  # count NAs per stations
  numNA <- apply(matrix(stfdf@data[,covar],
                        nrow=nrowsp,byrow= FALSE), MARGIN=1,
                 FUN=function(x) sum(is.na(x)))
  rem <- numNA != length(time)
  stfdf <-  stfdf[rem,drop= FALSE]
}

# Remove dates out of covariates
rm.days <- c()
for (t in 1:length(time)) {
  if(sum(complete.cases(stfdf[, t]@data)) == 0) {
    rm.days <- c(rm.days, t)
  }
}
if(!is.null(rm.days)){
  stfdf <- stfdf[,-rm.days]
}

### Example with STFDF and without parallel processing and without refitting of variogram
results <- cv.strk(data = stfdf,
                   obs.col = 1, # "tempc"
                   data.staid.x.y.z = c(1,NA,NA,NA), # first column is station ID
                   reg.coef = tregcoef[[1]],
                   vgm.model = tvgms[[1]],
                   sp.nmax = 20,
                   time.nmax = 2,
                   type = "LLO",
                   k = 5,
                   seed = 42,
                   refit = FALSE,
                   progress = TRUE
)

stplot(results[,,"pred"])
summary(results)
# accuracy
acc.metric.fun(results@data$obs, results@data$pred, "R2")
acc.metric.fun(results@data$obs, results@data$pred, "RMSE")
acc.metric.fun(results@data$obs, results@data$pred, "MAE")
acc.metric.fun(results@data$obs, results@data$pred, "CCC")

### Example with data.frame, parallel processing, and refit
stfdf.df <- as.data.frame(stfdf)
results <- cv.strk(data = stfdf.df,
                   obs.col = "tempc",
                   data.staid.x.y.z = c("staid","lon","lat","time"),
                   zero.tol = 0,
                   reg.coef = tregcoef[[1]],
                   vgm.model = tvgms[[1]],
                   sp.nmax = 20,
                   time.nmax = 2,
                   type = "LLO",
                   k = 5,
                   seed = 42,
                   refit = TRUE,
                   output.format = "STFDF",
                   parallel.processing = TRUE,
                   pp.type = "doParallel", # "snowfall"
                   cpus = detectCores()-1,
                   progress = TRUE
)

stplot(results[,,"pred"])
summary(results)
# accuracy
acc.metric.fun(results@data$obs, results@data$pred, "R2")
acc.metric.fun(results@data$obs, results@data$pred, "RMSE")
acc.metric.fun(results@data$obs, results@data$pred, "MAE")
acc.metric.fun(results@data$obs, results@data$pred, "CCC")

