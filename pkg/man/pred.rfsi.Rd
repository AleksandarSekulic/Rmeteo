\name{pred.rfsi}

\alias{pred.rfsi}

\title{Random Forest Spatial Interpolation (RFSI) prediction}

\description{Function for prediction based on Random Forest Spatial Interpolation (RFSI) model (Sekulić et al. 2020).}

\usage{
pred.rfsi(model,
          data,
          zcol=1,
          data.staid.x.y.time = c(1,2,3,4),
          obs,
          obs.staid.time = c(1,2),
          stations,
          stations.staid.x.y = c(1,2,3),
          newdata,
          newdata.staid.x.y.time = c(1,2,3),
          zero.tol = 0,
          s.crs = NA,
          newdata.s.crs=NA,
          t.crs = NA,
          output.format = "data.frame",
          cpus = detectCores()-1,
          progress = TRUE,
          ...)
}

\arguments{
  \item{model}{ranger; An RFSI model made by \link{rfsi} function.}
  \item{data}{\link[spacetime]{STFDF-class}, \link[spacetime]{STSDF-class} or \code{data.frame}; Contains dependent variable (observations) used for prediction from an RFSI model. If \code{data.frame} object, it should have next columns: station ID (staid), longitude (x), latitude (y), time of the observation (time), and observation value (obs).}
  \item{zcol}{numeric or character; Column name or number showing the position of dependent variable (observations) in \code{data}. Default is 1.}
  \item{data.staid.x.y.time}{numeric or character vector; Positions or names of the station ID (staid), longitude (x), latitude (y) and time columns in \code{data} if \code{data} is data.frame. Default is c(1,2,3,4).}
  \item{obs}{\code{data.frame}; Contains dependent variable (observations). It should have next columns: station ID (staid), time of the observation (time), and observation value (obs)). This object is used together with \code{stations} (see below) to create \link[spacetime]{STFDF-class} object (if \code{data} object is missing) which is then used for prediction from an RFSI model.}
  \item{obs.staid.time}{numeric or character vector; Positions or names of the station ID (staid) and time columns in \code{obs}. Default is c(1,2).}
  \item{stations}{\code{data.frame}; It should have next columns: station ID (staid), longitude (x) and latitude (y) of the stations. This object is used together with \code{obs} (see above) if \code{data} object is missing.}
  \item{stations.staid.x.y}{numeric or character vector; Positions or names of the station ID (staid), longitude (x) and latitude (y) columns in \code{stations}. Default is c(1,2,3).}
  \item{newdata}{\link[spacetime]{STFDF-class}, \link[spacetime]{STSDF-class} or \code{data.frame}; Contains prediction locations and covariates. If \code{data.frame} object, it should have next columns: prediction location ID (staid), longitude (x), latitude (y), time of the prediction (time) and covariates (cov1, cov2, ...). Covariate names should be the same as in the \code{model}.}
  \item{newdata.staid.x.y.time}{numeric or character vector; Positions or names of the prediction location ID (staid), longitude (x), latitude (y) and time columns in \code{newdata}. Default is c(1,2,3).}
  \item{zero.tol}{numeric; A distance value below (or equal to) which locations are considered as duplicates. Default is 0. See \link{rm.dupl}.}
  \item{s.crs}{Source \link[sp]{CRS} of observations (\code{data}). If NA, read from \code{data}}
  \item{newdata.s.crs}{Source \link[sp]{CRS} of (\code{newdata}). If NA, read from \code{data}}
  \item{t.crs}{Target \link[sp]{CRS} for observations (\code{data}) and (\code{newdata}) reprojection. Note that observations should be in projection for finding nearest observations based on Eucleadean distances (see function \link{near.obs}) and should be the same as one used for making an RFSI model.}
  \item{output.format}{character; Format of the output, "data.frame" (default), "SpatialPointsDataFrame", "SpatialPixelsDataFrame", "STFDF" or "STSDF".}
  \item{cpus}{numeric; Number of processing units. Default is detectCores()-1.}
  \item{progress}{logical; If progress bar is shown. Default is TRUE.}
  \item{...}{Further arguments passed to \link[ranger]{ranger}.}
}

\value{
A \code{data.frame}, \link[sp]{SpatialPointsDataFrame}, \link[sp]{SpatialPixelsDataFrame}, \link[spacetime]{STFDF-class} or \link[spacetime]{STSDF-class} object (depends on \code{output.format} argument) with prediction (\code{pred} column).
}

\author{Aleksandar Sekulic \email{asekulic@grf.bg.ac.rs}}

\references{
Sekulić, A., Kilibarda, M., Heuvelink, G. B., Nikolić, M. & Bajat, B. Random Forest Spatial Interpolation. Remote. Sens. 12, 1687, https://doi.org/10.3390/rs12101687 (2020).
}

\seealso{
  \code{\link{near.obs}}
  \code{\link{rfsi}}
  \code{\link{tune.rfsi}}
  \code{\link{cv.rfsi}}
}

\examples{
library(sp)
library(spacetime)
library(gstat)
library(plyr)
library(CAST)
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

rfsi_model <- rfsi(formula=formula,
                   data= stfdf,
                   zero.tol=0,
                   n.obs=5, # nearest obs
                   s.crs=NA, # read from data
                   t.crs=t.crs,
                   cpus=detectCores()-1,
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

# Make prediction from the RFSI model
newdata = regdata #[,1,drop=FALSE]
# newdata = regdata.df
# newdata.staid.x.y.time = c(1,2,3,4) # if data.frame
# spdf <- newdata[,1]
# spdf@data <- cbind(spdf@data, newdata@sp@data)
output.format = "STFDF" # STSDF # data.frame
zcol=1
model = rfsi_model
newdata.s.crs=NA

source('~/meteo/packageWV/pkg/R/pred.rfsi.R')
rfsi_prediction <- pred.rfsi(model = rfsi_model,
                             data = stfdf,
                             zcol=1,
                             newdata = newdata, # spdf,
                             # newdata.staid.x.y.time = newdata.staid.x.y.time,
                             output.format = "STFDF", # "SpatialPixelsDataFrame",
                             zero.tol=0,
                             s.crs=s.crs,
                             newdata.s.crs=newdata.s.crs,
                             t.crs=t.crs,
                             cpus=detectCores()-1,
                             progress=TRUE
)

stplot(rfsi_prediction)
# plot(rfsi_prediction['pred'])

}

