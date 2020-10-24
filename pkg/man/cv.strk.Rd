\name{cv.strk}

\alias{cv.strk}

\title{k-fold cross-validation for spatio-temporal regression kriging}

\description{k-fold cross-validation function for spatio-temporal regression kriging based on \link{pred.strk}. The cross-validation is made for \link[spacetime]{STFDF-class}, \link[spacetime]{STSDF-class}, or for data.frame objects.
Currently, only spatial (leave-location-out) cross-validation is implemented. Temporal and spatio-temporal cross-validation will be implemented in the future.}

\usage{
cv.strk(data,
        zcol=1,
        data.staid.x.y.time = c(1,2,3,4),
        obs,
        obs.staid.time = c(1,2),
        stations,
        stations.staid.x.y = c(1,2,3),
        zero.tol=0,
        reg.coef,
        vgm.model,
        sp.nmax=20,
        time.nmax=2,
        type = "LLO",
        k = 5,
        seed = 42,
        folds,
        fold.column,
        refit = TRUE,
        output.format = "STFDF",
        parallel.processing = FALSE,
        pp.type = "snowfall",
        cpus=detectCores()-1,
        progress=TRUE,
        ...)
}

\arguments{
  \item{data}{\link[spacetime]{STFDF-class}, \link[spacetime]{STSDF-class} or \code{data.frame}; Contains dependent variable (observations) and covariates in space and time used to perform cross-validation for spatio-temporal regression kriging. If \code{data.frame} object, it should have next columns: station ID (staid), longitude (x), latitude (y), time of the observation (time), observation value (obs) and covariates (cov1, cov2, ...). Covariate names should be the same as in the \code{reg.coef} (see below). If covariates are missing, then cross-validation for spatio-temporal ordinary kriging is performed.}
  \item{zcol}{numeric or character; Column name or number showing the position of dependent variable (observations) in \code{data}. Default is 1.}
  \item{data.staid.x.y.time}{numeric or character vector; Positions or names of the station ID (staid), longitude (x), latitude (y) and time columns in \code{data} if \code{data} is data.frame. Default is c(1,2,3,4).}
  \item{obs}{\code{data.frame}; Contains dependent variable (observations) and covariates in space and time. It should have next columns: station ID (staid), time of the observation (time), observation value (obs) and covariates (cov1, cov2, ...). This object is used together with \code{stations} (see below) to create \link[spacetime]{STFDF-class} object (if \code{data} object is missing) which is then used to perform cross-validation for spatio-temporal regression kriging.  If covariates are missing, cross-validation for spatio-temporal ordinary kriging is performed.}
  \item{obs.staid.time}{numeric or character vector; Positions or names of the station ID (staid) and time columns in \code{obs}. Default is c(1,2).}
  \item{stations}{\code{data.frame}; It should have next columns: station ID (staid), longitude (x) and latitude (y) of the stations. This object is used together with \code{obs} (see above) if \code{data} object is missing.}
  \item{stations.staid.x.y}{numeric or character vector; Positions or names of the station ID (staid), longitude (x) and latitude (y) columns in \code{stations}. Default is c(1,2,3).}
  \item{zero.tol}{numeric; A distance value below (or equal to) which locations are considered as duplicates. Default is 0. See \link{rm.dupl}. Duplicates are removed to avoid singular covariance matrices in kriging.}
 \item{reg.coef}{named numeric vector; Named linear regression coefficients. Names of the coefficients (e.g. "Intercept", "temp_geo", "modis", "dem", "twi") will be used to match appropriate covariates from \code{data} (or \code{obs}). Coefficients for metorological variables (temperature, precipitation, etc.) can be taken from data(\link{tregcoef}) or can be specified by the user.}
  \item{vgm.model}{StVariogramModel list; Spatio-temporal variogram of regression residuals (or observations if spatio-temporal ordinary kriging). See \link[gstat]{vgmST}. Spatio-temporal variogram model on residuals for metorological variables (temperature, precipitation, etc.) can be taken from data(\link{tvgms}) or can be specified by the user as a \link[gstat]{vgmST} object.}
  \item{sp.nmax}{numeric; A number of spatially nearest observations that should be used for fold kriging predictions. If /code{tiling} is TRUE (see below), then is a number of spatially nearest observations that should be used for each tile. Deafult is 20.}
  \item{time.nmax}{numeric; A number of temporally nearest observations that should be used for fold kriging predictions Deafult is 2.}
  \item{type}{character; Type of cross-validation: leave-location-out ("LLO"), leave-time-out ("LTO"), and leave-location-time-out ("LLTO"). Default is "LLO". "LTO" and "LLTO" are not implemented yet. Will be in the future.}
  \item{k}{numeric; Number of random folds that will be created with \link[CAST]{CreateSpacetimeFolds} function if \code{folds} or \code{fold.column} parameters are missing. Default is 5.}
  \item{seed}{numeric; Random seed that will be used to generate folds with \link[CAST]{CreateSpacetimeFolds} function.}
  \item{folds}{numeric or character vector; Showing folds of \code{data} observations. Used if \code{fold.column} parameter is not specified.}
  \item{fold.column}{numeric or character; Column name or number showing the position of variable in \code{data} that represents folds.}
  \item{refit}{logical; If refit of linear regression trend and spatio-teporal variogram should be performed. Spatio-teporal variogram is fit using \code{vgm.model} as desired spatio-temporal model for \link[gstat]{fit.StVariogram} function. Default is TRUE.}
  \item{output.format}{character; Format of the output, "STFDF" (default), "STSDF" or "data.frame" (\code{data.frame}).}
  \item{parallel.processing}{logical; If parallel processing is performed. Default is FALSE.}
  \item{pp.type}{character; Type (R package) of parallel processing, "snowfall" (default) or "doParallel".}
  \item{cpus}{numeric; Number of processing units. Default is detectCores()-1.}
  \item{progress}{logical; If progress bar is shown. Default is TRUE.}
  \item{...}{Further arguments passed to \link[gstat]{krigeST} or \link{pred.strk}.}
}

\value{
A \link[spacetime]{STFDF-class}, \link[spacetime]{STSDF-class} or \code{data.frame} obejct (depends on \code{output.format} argument), with columns:
  \item{obs}{Observations.}
  \item{pred}{Predictions from cross-validation.}
  \item{folds}{Folds used for cross-validation.}

For accuracy metrics see \link{acc.metric.fun} function.
}

\author{Aleksandar Sekulic \email{asekulic@grf.bg.ac.rs}, Milan Kilibarda \email{kili@grf.bg.ac.rs}}

\references{
Kilibarda, M., T. Hengl, G. B. M. Heuvelink, B. Graeler, E. Pebesma, M. Percec Tadic, and B. Bajat (2014), Spatio-temporal interpolation of daily temperatures for global land areas at 1 km resolution, J. Geophys. Res. Atmos., 119, 2294-2313, doi:10.1002/2013JD020803.
}

\seealso{
 \code{\link{acc.metric.fun}}
 \code{\link{pred.strk}}
 \code{\link{tregcoef}}
 \code{\link{tvgms}}
 \code{\link{regdata}}
 \code{\link{meteo2STFDF}}
 \code{\link{tgeom2STFDF}}
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

# Cross-validation example with data.frames, parallel processing, and refit
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
# accuracy
acc.metric.fun(results@data$obs, results@data$pred, "RMSE")

}

