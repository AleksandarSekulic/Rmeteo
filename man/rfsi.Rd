\name{rfsi}

\alias{rfsi}

\title{Random Forest Spatial Interpolation (RFSI) model}

\description{Function for creation of Random Forest Spatial Interpolation (RFSI) model (Sekulić et al. 2020). Besides environmental covariates, RFSI uses additional spatial covariates: (1) observations at n nearest locations and (2) distances to them, in order to include spatial context into the random forest.}

\usage{
rfsi(formula,
     data,
     data.staid.x.y.time = c(1,2,3,4),
     obs,
     obs.staid.time = c(1,2),
     stations,
     stations.staid.x.y = c(1,2,3),
     zero.tol = 0,
     n.obs = 10,
     avg = FALSE,
     increment,
     range,
     direct = FALSE,
     use.idw = FALSE,
     idw.p = 2,
     s.crs = NA,
     t.crs = NA,
     cpus = detectCores()-1,
     progress = TRUE,
     ...)
}

\arguments{
  \item{formula}{formula; Specifying the dependent variable (without nearest observations and distances to them). If \code{z~1}, an RFSI model will be made from nearest obsevrations and distances to them as covariates.}
  \item{data}{\link[spacetime]{STFDF-class}, \link[spacetime]{STSDF-class} or \code{data.frame}; Contains dependent variable (observations) and covariates used for making an RFSI model. If \code{data.frame} object, it should have next columns: station ID (staid), longitude (x), latitude (y), time of the observation (time), observation value (obs) and covariates (cov1, cov2, ...). If covariates are missing, the RFSI model will be made from nearest obsevrations and distances to them as covariates (\code{formula=z~1}).}
  \item{data.staid.x.y.time}{numeric or character vector; Positions or names of the station ID (staid), longitude (x), latitude (y) and time columns in \code{data} if \code{data} is data.frame. Default is c(1,2,3,4).}
  \item{obs}{\code{data.frame}; Contains dependent variable (observations) and covariates in space and time. It should have next columns: station ID (staid), time of the observation (time), observation value (obs) and covariates (cov1, cov2, ...). This object is used together with \code{stations} (see below) to create \link[spacetime]{STFDF-class} object (if \code{data} object is missing) which is then used for making an RFSI model.  If covariates are missing, the RFSI model will be made from nearest obsevrations and distances to them as covariates (\code{formula=z~1}).}
  \item{obs.staid.time}{numeric or character vector; Positions or names of the station ID (staid) and time columns in \code{obs}. Default is c(1,2).}
  \item{stations}{\code{data.frame}; It should have next columns: station ID (staid), longitude (x) and latitude (y) of the stations. This object is used together with \code{obs} (see above) if \code{data} object is missing.}
  \item{stations.staid.x.y}{numeric or character vector; Positions or names of the station ID (staid), longitude (x) and latitude (y) columns in \code{stations}. Default is c(1,2,3).}
  \item{zero.tol}{numeric; A distance value below (or equal to) which locations are considered as duplicates. Default is 0. See \link{rm.dupl}.}
  \item{n.obs}{numeric; Number of nearest observations to be used as covariates in RFSI model (see function \link{near.obs}). Note that it cannot be larger than number of obsevrations. Default is 10.}
  \item{avg}{boolean; Will averages in circles with different radiuses be calculated and used as covariates (see function \link{near.obs}). Default is FALSE.}
  \item{increment}{numeric; Increment of radiuses for calculation of averages in circles with different radiuses (see function \link{near.obs}). Units depends on CRS of coordinates.}
  \item{range}{numeric; Maximum radius for calculation of averages in circles with different radiuses (see function \link{near.obs}). Units depends on CRS of coordinates.}
  \item{direct}{boolean; Will nearest observation in quadrants be calculated and used as covariates (see function \link{near.obs}). Default is FALSE.}
  \item{use.idw}{boolean; Will IDW predictions from \code{n.obs} nearest observations be calculated and used as covariate (see function \link{near.obs}). Default is FALSE.}
  \item{idw.p}{numeric; Exponent for IDW weights (see function \link{near.obs}). Default is 2.}
  \item{s.crs}{Source \link[sp]{CRS-class} of observations (\code{data}). If NA, read from \code{data}}
  \item{t.crs}{Target \link[sp]{CRS-class} for observations (\code{data}) reprojection. If NA, will be set to \code{s.crs}. Note that observations should be in projection for finding nearest observations based on Eucleadean distances (see function \link{near.obs}).}
  \item{cpus}{numeric; Number of processing units. Default is detectCores()-1.}
  \item{progress}{logical; If progress bar is shown. Default is TRUE.}
  \item{...}{Further arguments passed to \link[ranger]{ranger}.}
}

\value{
An RFSI model of class \link[ranger]{ranger}.
}

\author{Aleksandar Sekulic \email{asekulic@grf.bg.ac.rs}}

\references{
Sekulić, A., Kilibarda, M., Heuvelink, G. B., Nikolić, M. & Bajat, B. Random Forest Spatial Interpolation. Remote. Sens. 12, 1687, https://doi.org/10.3390/rs12101687 (2020).
}

\seealso{
  \code{\link{near.obs}}
  \code{\link{pred.rfsi}}
  \code{\link{tune.rfsi}}
  \code{\link{cv.rfsi}}
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

# Make an RFSI model
formula = 'tempc ~ temp_geo + modis + dem + twi'  # without nearest obs
t.crs=CRS("+proj=utm +zone=34 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
s.crs=NA

# ### obs + stations
# obs = as.data.frame(stfdf) # data.frame(id,time,obs,cov)
# obs.staid.time = c(3,4)
# stations = as.data.frame(stfdf) # data.frame(id,x,y)
# stations.staid.x.y = c(3,1,2)
# s.crs=CRS("+proj=longlat +datum=WGS84")
# 
# ### data.frame
# data1 = as.data.frame(stfdf)
# data.staid.x.y.time <- c(3,1,2,4)
# s.crs=CRS("+proj=longlat +datum=WGS84")
# 
# ### spatial only
# # SPDF
# data1 = as.data.frame(stfdf[,1])
# data1$staid = 1:10
# data.staid.x.y.time <- c(8,6,7,NA)
# 
# obs = data1 # data.frame(id,time,obs,cov)
# obs.staid.time = c(8,NA)
# stations = data1 # data.frame(id,x,y)
# stations.staid.x.y = c(8,6,7)
###

rfsi_model <- rfsi(formula=formula,
                data= stfdf, # data1,
                # data.staid.x.y.time = data.staid.x.y.time,
                # obs = obs,
                # obs.staid.time = obs.staid.time,
                # stations = stations,
                # stations.staid.x.y = stations.staid.x.y,
                zero.tol=0,
                n.obs=5, # nearest obs
                s.crs=s.crs,
                t.crs=t.crs,
                cpus=2, # detectCores()-1,
                progress=TRUE,
                # ranger parameters
                importance = "impurity",
                seed = 42,
                num.trees = 250,
                mtry = 5,
                splitrule = "variance",
                min.node.size = 5,
                sample.fraction = 0.95,
                quantreg = FALSE)

rfsi_model
# Type:                             Regression 
# Number of trees:                  250 
# Sample size:                      20 
# Number of independent variables:  14 
# Mtry:                             5 
# Target node size:                 5 
# Variable importance mode:         impurity 
# Splitrule:                        variance 
# OOB prediction error (MSE):       0.5238891 
# R squared (OOB):                  0.3167993  
sort(rfsi_model$variable.importance)
sum("obs" == substr(rfsi_model$forest$independent.variable.names, 1, 3))

}

