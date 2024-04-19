\name{cv.strk}

\alias{cv.strk}

\title{k-fold cross-validation for spatio-temporal regression kriging}

\description{k-fold cross-validation function for spatio-temporal regression kriging based on \link{pred.strk}. Currently, only spatial (leave-location-out) cross-validation is implemented. Temporal and spatio-temporal cross-validation will be implemented in the future.}

\usage{
cv.strk(data,
        obs.col=1,
        data.staid.x.y.z = NULL,
        crs = NA,
        zero.tol=0,
        reg.coef,
        vgm.model,
        sp.nmax=20,
        time.nmax=2,
        type = "LLO",
        k = 5,
        seed = 42,
        folds,
        refit = TRUE,
        output.format = "STFDF",
        parallel.processing = FALSE,
        pp.type = "snowfall",
        cpus=detectCores()-1,
        progress=TRUE,
        ...)
}

\arguments{
  \item{data}{\link[spacetime]{STFDF-class}, \link[spacetime]{STSDF-class}, \link[spacetime]{STIDF-class}, \link[sf:sf]{sf-class}, \link[sftime:st_sftime]{sftime-class}, \link[terra]{SpatVector-class} or \link[base]{data.frame}; Contains target variable (observations) and covariates in space and time used to perform STRK cross validation. If \link[base]{data.frame} object, it should have next columns: station ID (staid), longitude (x), latitude (y), 3rd component - time, depth, ... (z) of the observation, observation value (obs), and covariates (cov1, cov2, ...). Covariate names should be the same as in the \code{reg.coef} (see below). If covariates are missing, then spatio-temporal ordinary kriging cross validation is performed.}
  \item{obs.col}{numeric or character; Column name or number showing position of the observation column in the \code{data}. Default is 1.}
  \item{data.staid.x.y.z}{numeric or character vector; Positions or names of the station ID (staid), longitude (x), latitude (y) and 3rd component - time, depth (z) columns in \link[base]{data.frame} object (e.g. c(1,2,3,4)). If \code{data} is \link[sf:sf]{sf-class}, \link[sftime:st_sftime]{sftime-class}, or \link[terra]{SpatVector-class} object, \code{data.staid.x.y.z} is used to point staid and z position. If \code{data} is \link[spacetime]{STFDF-class}, \link[spacetime]{STSDF-class}, \link[spacetime]{STIDF-class} object, \code{data.staid.x.y.z} is used to point only staid position. Default is NULL.}
  \item{crs}{\link[sf]{st_crs} or \link[terra]{crs}; Source CRS of \code{data}. If \code{data} contains crs, \code{crs} will not be used. Default is NA.}
  \item{zero.tol}{numeric; A distance value below (or equal to) which locations are considered as duplicates. Default is 0. See \link{rm.dupl}. Duplicates are removed to avoid singular covariance matrices in kriging.}
  \item{reg.coef}{numeric; Vector of named linear regression coefficients. Names of the coefficients (e.g. "Intercept", "temp_geo", "modis", "dem", "twi") will be used to match appropriate covariates from \code{data}. Coefficients for metorological variables (temperature, precipitation, etc.) can be taken from data(\link{tregcoef}) or can be specified by the user.}
  \item{vgm.model}{StVariogramModel list; Spatio-temporal variogram of regression residuals (or observations if spatio-temporal ordinary kriging). See \link[gstat]{vgmST}. Spatio-temporal variogram model on residuals for metorological variables (temperature, precipitation, etc.) can be taken from data(\link{tvgms}) or can be specified by the user as a \link[gstat]{vgmST} object.}
  \item{sp.nmax}{numeric; A number of spatially nearest observations that should be used for kriging predictions. If \code{tiling} is TRUE (see below), then is a number of spatially nearest observations that should be used for each tile. Deafult is 20.}
  \item{time.nmax}{numeric; A number of temporally nearest observations that should be used for kriging predictions Deafult is 2.}
  \item{type}{character; Type of cross-validation: leave-location-out ("LLO"), leave-time-out ("LTO"), and leave-location-time-out ("LLTO"). Default is "LLO". "LTO" and "LLTO" are not implemented yet. Will be in the future.}
  \item{k}{numeric; Number of random folds that will be created with \link[CAST]{CreateSpacetimeFolds} function. Default is 5.}
  \item{seed}{numeric; Random seed that will be used to generate outer  and inner folds with \link[CAST]{CreateSpacetimeFolds} function.}
  \item{folds}{numeric or character vector or value; Showing folds column (if value) or rows (vector) of \code{data} observations used for cross-validation. If missing, will be created with \link[CAST]{CreateSpacetimeFolds} function.}
  \item{refit}{logical; If refit of linear regression trend and spatio-teporal variogram should be performed. Spatio-teporal variogram is fit using \code{vgm.model} as desired spatio-temporal model for \link[gstat]{fit.StVariogram} function. Default is TRUE.}
  \item{output.format}{character; Format of the output, \link[spacetime]{STFDF-class} (default), \link[spacetime]{STSDF-class}, \link[spacetime]{STIDF-class}, \link[base]{data.frame}, \link[sf:sf]{sf-class}, \link[sftime:st_sftime]{sftime-class}, or \link[terra]{SpatVector-class}.}
  \item{parallel.processing}{logical; If parallel processing is performed. Default is FALSE.}
  \item{pp.type}{character; Type (R package) of parallel processing, "snowfall" (default) or "doParallel".}
  \item{cpus}{numeric; Number of processing units. Default is detectCores()-1.}
  \item{progress}{logical; If progress bar is shown. Default is TRUE.}
  \item{...}{Further arguments passed to \link[gstat]{krigeST} or \link{pred.strk}.}
}

\value{
A \link[spacetime]{STFDF-class} (default), \link[spacetime]{STSDF-class}, \link[spacetime]{STIDF-class}, \link[base]{data.frame}, \link[sf:sf]{sf-class}, \link[sftime:st_sftime]{sftime-class}, or \link[terra]{SpatVector-class} object (depends on \code{output.format} argument), with columns:
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
                   data.staid.x.y.z = c(1,NA,NA,NA),
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

\donttest{
stplot(results[,,"pred"])
}
summary(results)
# accuracy
acc.metric.fun(results@data$obs, results@data$pred, "R2")
acc.metric.fun(results@data$obs, results@data$pred, "RMSE")
acc.metric.fun(results@data$obs, results@data$pred, "MAE")
acc.metric.fun(results@data$obs, results@data$pred, "CCC")

}

