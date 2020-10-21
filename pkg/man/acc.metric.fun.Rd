\name{acc.metric.fun}

\alias{acc.metric.fun}

\title{Accuracy metrics calculation}

\description{A function for calculation of classification and regression accuracy metrics.}

\usage{
acc.metric.fun(obs, pred, acc.m)
}

\arguments{
  \item{obs}{\code{numeric or factor vector;} Observations.}
  \item{pred}{\code{numeric or factor vector;} Predictions.}
  \item{acc.m}{\code{character}; Accuracy metric. Possible values for regression: "ME", "MAE", "RMSE", "R2", "CCC". Possible values for classification: "Accuracy", "Kappa", "AccuracyLower", "AccuracyUpper", "AccuracyNull", "AccuracyPValue", "McnemarPValue".}
}

\value{
Accuracy metric value.
}

\author{Aleksandar Sekulic \email{asekulic@grf.bg.ac.rs}}

\references{
Sekulić, A., Kilibarda, M., Heuvelink, G. B., Nikolić, M. & Bajat, B. Random Forest Spatial Interpolation.Remote. Sens. 12, 1687, https://doi.org/10.3390/rs12101687 (2020).
}

\seealso{
 \code{\link{acc.metric.fun}}
 \code{\link{rfsi}}
 \code{\link{pred.rfsi}}
 \code{\link{tune.rfsi}}
 \code{\link{cv.rfsi}}
 \code{\link{pred.strk}}
 \code{\link{cv.strk}}
}

\examples{
library(sp)
library(spacetime)
library(gstat)
library(plyr)
library(CAST)
library(doParallel)
library(ranger)
# prepare data
# load observation - data.frame of mean temperatures
data(dtempc)
data(stations)

serbia= point.in.polygon(stations$lon, stations$lat, c(18,22.5,22.5,18), c(40,40,46,46))
st= stations[ serbia!=0, ]
dtempc <- dtempc[dtempc$staid \%in\% st$staid, ]
dtempc <- dtempc[complete.cases(dtempc),]

# create STFDF
stfdf <- meteo2STFDF(dtempc,st)
# Adding CRS
stfdf@sp@proj4string <- CRS('+proj=longlat +datum=WGS84')

# load covariates for mean temperatures
data(regdata)
# str(regdata)
regdata@sp@proj4string <- CRS('+proj=longlat +datum=WGS84')

# load precalculated variograms
data(tvgms)
data(tregcoef)

# Overlay observations with covariates
time <- index(stfdf@time)
covariates.df <- as.data.frame(regdata)
names_covar <- names(tregcoef[[1]])[-1]
for (covar in names_covar){
  nrowsp <- length(stfdf@sp)
  regdata@sp=as(regdata@sp,'SpatialPixelsDataFrame')
  ov <- sapply(time, function(i) 
    if (covar \%in\% names(regdata@data)) {
      if (as.Date(i) \%in\% as.Date(index(regdata@time))) {
        over(stfdf@sp, as(regdata[, i, covar], 'SpatialPixelsDataFrame'))[, covar]
      } else {
        rep(NA, length(stfdf@sp))
      }
    } else {
      over(stfdf@sp, as(regdata@sp[covar], 'SpatialPixelsDataFrame'))[, covar]
    }
  )
  # ov <- do.call('cbind', ov)
  ov <- as.vector(ov)
  if (all(is.na(ov))) {
    stop(paste('There is no overlay of data with covariates!', sep = ""))
  }
  stfdf@data[covar] <- ov
}

# remove stations out of covariates
for (covar in names_covar){
  # count NAs per stations
  numNA <- apply(matrix(stfdf@data[,covar],
                        nrow=nrowsp,byrow=FALSE), MARGIN=1,
                 FUN=function(x) sum(is.na(x)))
  # Remove stations out of covariates
  rem <- numNA != length(time)
  stfdf <-  stfdf[rem,drop=FALSE]
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

# Cross-validation for mean temperature for days "2011-07-05" and "2011-07-06" 
# global model is used for regression and variogram

results <- cv.strk(data = stfdf,
                   zcol=1, # "tempc"
                   reg.coef = tregcoef[[1]],
                   vgm.model = tvgms[[1]],
                   sp.nmax = 20,
                   time.nmax = 2,
                   type = "LLO",
                   k = 5,
                   seed = 42,
                   refit = FALSE,
                   progress=TRUE
)

stplot(results)
summary(results)
# accuracy
acc.metric.fun(results@data$obs, results@data$pred, "RMSE")

# Cross-validatiin example with data.frames, parallel processing, and refit
# global model is used for regression, variogram, and refit

library(snowfall)
library(doParallel)
# create data.frames
stfdf.df <- as.data.frame(stfdf)

results <- cv.strk(data = stfdf.df,
                   zcol="tempc",
                   data.staid.x.y.time = c("staid","lon","lat","time"),
                   # obs = stfdf.df, # if used, comment data argument
                   # obs.staid.time = c("staid","time"),
                   # stations = stfdf.df,
                   # stations.staid.x.y = c("staid","lon","lat"),
                   zero.tol=0,
                   reg.coef = tregcoef[[1]],
                   vgm.model = tvgms[[1]],
                   sp.nmax = 20,
                   time.nmax = 2,
                   type = "LLO",
                   k = 5,
                   seed = 42,
                   refit = TRUE,
                   parallel.processing = TRUE,
                   pp.type = "snowfall", # "doParallel"
                   cpus=detectCores()-1,
                   progress=TRUE
)

stplot(results)
summary(results)
# accuracy metric calculation
acc.metric.fun(results@data$obs, results@data$pred, "RMSE")

}

