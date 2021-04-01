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
data(dtempc) # temperature data
data(stations) # station locations
data(regdata) # covariates
regdata@sp@proj4string <- CRS('+proj=longlat +datum=WGS84')
data(tvgms) # ST variogram models
data(tregcoef) # MLR coefficients
regdata.df <- as.data.frame(regdata)

serbia= point.in.polygon(stations$lon, stations$lat, c(18,22.5,22.5,18), c(40,40,46,46))
st= stations[ serbia!=0, ]
dtempc <- dtempc[dtempc$staid %in% st$staid, ]
dtempc <- dtempc[complete.cases(dtempc),]
summary(dtempc)
# create data.frame
stfdf.df <- join(dtempc, st)
summary(stfdf.df)
# create STFDF
stfdf <- meteo2STFDF(dtempc,st)
stfdf@sp@proj4string <- CRS('+proj=longlat +datum=WGS84')

#################### pred.strk ####################

# Calculate prediction of mean temperatures for "2011-07-05" and "2011-07-06" 
# global model is used for regression and variogram

### Example with STFDF and without parallel processing
results <- pred.strk(data = stfdf, # observations
                     newdata = regdata, # prediction locations with covariates
                     # newdata = regdata[,2,drop=FALSE], # for one day only
                     output.format = "stfdf",
                     reg.coef = tregcoef[[1]], # MLR coefficients
                     vgm.model = tvgms[[1]], # STRK variogram model
                     sp.nmax = 20,
                     time.nmax = 2,
                     computeVar=TRUE
)

# plot prediction
stplot(results[,,"pred", drop=F], col.regions=bpy.colors())
stplot(results[,,"var", drop=F], col.regions=bpy.colors())

### Example with data.frames and parallel processing
results <- pred.strk(data = stfdf.df,
                     zcol = 3,
                     data.staid.x.y.time = c(1,4,5,2),
                     # obs = stfdf.df, # if used, comment data argument
                     # obs.staid.time = c(1,2),
                     # stations = stfdf.df,
                     # stations.staid.x.y = c(1,4,5),
                     newdata = regdata.df,
                     # newdata = regdata.df[regdata.df$time=="2011-07-06", ], # for one day only
                     newdata.staid.x.y.time = c(3,1,2,4),
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
stplot(results[,,"pred", drop=F], col.regions=bpy.colors())
stplot(results[,,"var", drop=F], col.regions=bpy.colors())

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
                        nrow=nrowsp,byrow=F), MARGIN=1,
                 FUN=function(x) sum(is.na(x)))
  rem <- numNA != length(time)
  stfdf <-  stfdf[rem,drop=F]
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
                   zcol = 1, # "tempc"
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
                   zcol = "tempc",
                   data.staid.x.y.time = c("staid","lon","lat","time"),
                   # obs = stfdf.df, # if used, comment data argument
                   # obs.staid.time = c("staid","time"),
                   # stations = stfdf.df,
                   # stations.staid.x.y = c("staid","lon","lat"),
                   zero.tol = 0,
                   reg.coef = tregcoef[[1]],
                   vgm.model = tvgms[[1]],
                   sp.nmax = 20,
                   time.nmax = 2,
                   type = "LLO",
                   k = 5,
                   seed = 42,
                   refit = TRUE,
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
