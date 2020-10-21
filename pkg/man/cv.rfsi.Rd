\name{cv.rfsi}

\alias{cv.rfsi}

\title{Nested k-fold cross-validation for Random Forest Spatial Interpolation (RFSI)}

\description{Function for nested k-fold cross-validation function for Random Forest Spatial Interpolation (RFSI) (Sekulić et al. 2020). It is based on \link{rfsi}, \link{pred.rfsi}, and \link{tune.rfsi} functions. Currently, only spatial (leave-location-out) cross-validation is implemented. Temporal and spatio-temporal cross-validation will be implemented in the future.}

\usage{
cv.rfsi(formula,
        data,
        data.staid.x.y.time = c(1,2,3,4),
        obs,
        obs.staid.time = c(1,2),
        stations,
        stations.staid.x.y = c(1,2,3),
        zero.tol = 0,
        use.idw = FALSE,
        s.crs = NA,
        t.crs = NA,
        tgrid,
        tgrid.n=10,
        tune.type = "LLO",
        k = 5,
        seed=42,
        folds,
        fold.column,
        acc.metric,
        output.format = "data.frame",
        cpus = detectCores()-1,
        progress = TRUE,
        ...)
}

\arguments{
  \item{formula}{formula; Specifying the dependent variable (without nearest observations and distances to them). If \code{z~1}, an RFSI model with nearest obsevrations and distances to them as covariates will be cross-validated.}
  \item{data}{\link[spacetime]{STFDF-class}, \link[spacetime]{STSDF-class} or \code{data.frame}; Contains dependent variable (observations) and covariates used for making an RFSI model. If \code{data.frame} object, it should have next columns: station ID (staid), longitude (x), latitude (y), time of the observation (time), observation value (obs) and covariates (cov1, cov2, ...). If covariates are missing, the RFSI model with nearest obsevrations and distances to them as covariates (\code{formula=z~1}) will be cross-validated.}
  \item{data.staid.x.y.time}{numeric or character vector; Positions or names of the station ID (staid), longitude (x), latitude (y) and time columns in \code{data} if \code{data} is data.frame. Default is c(1,2,3,4).}
  \item{obs}{\code{data.frame}; Contains dependent variable (observations) and covariates in space and time. It should have next columns: station ID (staid), time of the observation (time), observation value (obs) and covariates (cov1, cov2, ...). This object is used together with \code{stations} (see below) to create \link[spacetime]{STFDF-class} object (if \code{data} object is missing) which is then used for making an RFSI model.  If covariates are missing, the RFSI model with nearest obsevrations and distances to them as covariates (\code{formula=z~1}) will be cross-validated.}
  \item{obs.staid.time}{numeric or character vector; Positions or names of the station ID (staid) and time columns in \code{obs}. Default is c(1,2).}
  \item{stations}{\code{data.frame}; It should have next columns: station ID (staid), longitude (x) and latitude (y) of the stations. This object is used together with \code{obs} (see above) if \code{data} object is missing.}
  \item{stations.staid.x.y}{numeric or character vector; Positions or names of the station ID (staid), longitude (x) and latitude (y) columns in \code{stations}. Default is c(1,2,3).}
  \item{zero.tol}{numeric; A distance value below (or equal to) which locations are considered as duplicates. Default is 0. See \link{rm.dupl}.}
  \item{use.idw}{boolean; Will IDW predictions from \code{n.obs} nearest observations be calculated and used as covariate (see function \link{near.obs}). Default is FALSE.}
  \item{s.crs}{Source \link[sp]{CRS} of observations (\code{data}). If NA, read from \code{data}}
  \item{t.crs}{Target \link[sp]{CRS} for observations (\code{data}) reprojection. If NA, will be set to \code{s.crs}. Note that observations should be in projection for finding nearest observations based on Eucleadean distances (see function \link{near.obs}).}
  \item{tgrid}{data.frame; Possible tuning parameters for nested folds. The columns are named the same as the tuning parameters. Possible tuning parameters are: \code{n.obs}, \code{num.trees}, \code{mtry}, \code{min.node.size}, \code{sample.fraction}, \code{splirule}, and \code{idw.p}.}
  \item{tgrid.n}{numeric; number of randomly chosen \code{tgrid} combinations used for nested tuning of RFSI. If larger than \code{tgrid}, will be set to \code{length(tgrid)}}
  \item{tune.type}{character; Type of cross-validation used for nested tuning: leave-location-out ("LLO"), leave-time-out ("LTO"), and leave-location-time-out ("LLTO"). Default is "LLO". "LTO" and "LLTO" are not implemented yet. Will be in the future.}
  \item{k}{numeric; Number of random folds for cross-validation and nested tuning that will be created, with \link[CAST]{CreateSpacetimeFolds} function if \code{folds} or \code{fold.column} parameters are missing. Default is 5.}
  \item{seed}{numeric; Random seed that will be used to generate folds for cross-validation and nested tuning, with \link[CAST]{CreateSpacetimeFolds} function.}
  \item{folds}{numeric or character vector; Showing folds of \code{data} observations used for cross-validation. Used if \code{fold.column} parameter is not specified. Note that folds for nested tuning will be created with \link[CAST]{CreateSpacetimeFolds} function.}
  \item{fold.column}{numeric or character; Column name or number showing the position of variable in \code{data} that represents foldsused for cross-validation. Note that folds for nested tuning will be created with \link[CAST]{CreateSpacetimeFolds} function.}
  \item{acc.metric}{character; Accuracy metric that will be used as a criteria for choosing an optimal RFSI model in nested tuning. Possible values for regression: "ME", "MAE", "RMSE" (default), "R2", "CCC". Possible values for classification: "Accuracy","Kappa" (default), "AccuracyLower", "AccuracyUpper", "AccuracyNull", "AccuracyPValue", "McnemarPValue".}
  \item{output.format}{character; Format of the output, "STFDF" (default), "STSDF" or "data.frame" (\code{data.frame}).}
  \item{cpus}{numeric; Number of processing units. Default is detectCores()-1.}
  \item{progress}{logical; If progress bar is shown. Default is TRUE.}
  \item{...}{Further arguments passed to \link[ranger]{ranger}.}
}

\value{
A \link[spacetime]{STFDF-class}, \link[spacetime]{STSDF-class} or \code{data.frame} obejct (depends on \code{output.format} argument), with columns:
  \item{obs}{Observations.}
  \item{pred}{Predictions from cross-validation.}
  \item{folds}{Folds used for cross-validation.}
}

\author{Aleksandar Sekulic \email{asekulic@grf.bg.ac.rs}}

\references{
Sekulić, A., Kilibarda, M., Heuvelink, G. B., Nikolić, M. & Bajat, B. Random Forest Spatial Interpolation. Remote. Sens. 12, 1687, https://doi.org/10.3390/rs12101687 (2020).
}

\seealso{
  \code{\link{near.obs}}
  \code{\link{rfsi}}
  \code{\link{pred.rfsi}}
  \code{\link{tune.rfsi}}
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
data(tregcoef)
# str(regdata)
regdata@sp@proj4string <- CRS('+proj=longlat +datum=WGS84')

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

formula = 'tempc ~ temp_geo + modis + dem + twi'  # without nearest obs
t.crs=CRS("+proj=utm +zone=34 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
s.crs=NA

# making tgrid
n.obs <- 2:3
min.node.size <- 2:10
sample.fraction <- seq(1, 0.632, -0.05) # 0.632 without / 1 with replacement
splitrule <- "variance"
ntree <- 250 # 500
mtry <- 3:(2+2*max(n.obs))
tgrid = expand.grid(min.node.size=min.node.size, num.trees=ntree,
                    mtry=mtry, n.obs=n.obs, sample.fraction=sample.fraction)

# Cross-validation of RFSI
rfsi_cv <- cv.rfsi(formula=formula, # without nearest obs
                   data= stfdf,
                   zero.tol=0,
                   s.crs=s.crs,
                   t.crs=t.crs,
                   tgrid=tgrid,
                   tgrid.n=5,
                   tune.type = "LLO",
                   k = 3, # number of folds
                   seed = 42,
                   acc.metric = "RMSE",
                   output.format = "data.frame", # "STFDF",
                   cpus=detectCores()-1,
                   progress=FALSE,
                   # ranger parameters
                   importance = "impurity")

summary(rfsi_cv)
# stplot(rfsi_cv[, , "pred"])
# stplot(rfsi_cv[, , "obs"])

}

