\name{pred.rfsi}

\alias{pred.rfsi}

\title{Random Forest Spatial Interpolation (RFSI) prediction}

\description{Function for spatial/spatio-temporal prediction based on Random Forest Spatial Interpolation (RFSI) model (Sekulić et al. 2020).}

\usage{
pred.rfsi(model,
          data,
          obs.col=1,
          data.staid.x.y.z = NULL,
          newdata,
          newdata.staid.x.y.z = NULL,
          zero.tol = 0,
          s.crs = NA,
          newdata.s.crs=NA,
          p.crs = NA,
          output.format = "data.frame",
          cpus = detectCores()-1,
          progress = TRUE,
          soil3d = FALSE, # soil RFSI
          depth.range = 0.1, # in units of depth
          no.obs = 'increase',
          ...)
}

\arguments{
  \item{model}{ranger; An RFSI model made by \link{rfsi} function.}
  \item{data}{\link[sf]{sf-class}, \link[sftime]{sftime-class}, \link[terra]{SpatVector-class} or \code{data.frame}; Contains target variable (observations) and covariates used for RFSI prediction. If \code{data.frame} object, it should have next columns: station ID (staid), longitude (x), latitude (y), 3rd component - time, depth, ... (z) of the observation, and observation value (obs).}
  \item{obs.col}{numeric or character; Column name or number showing position of the observation column in the \code{data}. Default is 1.}
  \item{data.staid.x.y.z}{numeric or character vector; Positions or names of the station ID (staid), longitude (x), latitude (y) and 3rd component (z) columns in data.frame \code{data} object (e.g. c(1,2,3,4)). If \code{data} is \link[sf]{sf-class}, \link[sftime]{sftime-class}, or \link[terra]{SpatVector-class} object, \code{data.staid.x.y.z} is used to point staid and z position. Set z position to NA (e.g. c(1,2,3,NA)) or ommit it (e.g. c(1,2,3)) for spatial interpolation. Default is NULL.}
  \item{newdata}{\link[sf]{sf-class}, \link[sftime]{sftime-class}, \link[terra]{SpatVector-class}, \link[terra]{SpatRaster-class} or \code{data.frame}; Contains prediction locations and covariates used for RFSI prediction. If \code{data.frame} object, it should have next columns: prediction location ID (staid), longitude (x), latitude (y), 3rd component - time, depth, ... (z), and covariates (cov1, cov2, ...). Covariate names have to be the same as in the \code{model}.}
  \item{newdata.staid.x.y.z}{numeric or character vector; Positions or names of the prediction location ID (staid), longitude (x), latitude (y) and 3rd component (z) columns in data.frame \code{newdata} object (e.g. c(1,2,3,4)). If \code{newdata} is \link[sf]{sf-class}, \link[sftime]{sftime-class}, \link[terra]{SpatVector-class} or \link[terra]{SpatRaster-class} object, \code{newdata.staid.x.y.z} is used to point staid and z position. Set z position to NA (e.g. c(1,2,3,NA)) or ommit it (e.g. c(1,2,3)) for spatial interpolation. Default is NULL.}
  \item{zero.tol}{numeric; A distance value below (or equal to) which locations are considered as duplicates. Default is 0. See \link{rm.dupl}.}
  \item{s.crs}{\link[sf]{crs-class}; Source CRS of \code{data}. If \code{data} contains crs, \code{s.crs} will be overwritten. Default is NA.}
  \item{newdata.s.crs}{\link[sf]{crs-class}; Source CRS of \code{newdata}. If \code{newdata} contains crs, \code{newdata.s.crs} will be overwritten. Default is NA.}
  \item{p.crs}{\link[sf]{crs-class}; Projection CRS for \code{data} reprojection. If NA, \code{s.crs} will be used for distance calculation. Note that observations should be in projection for finding nearest observations based on Eucleadean distances (see function \link{near.obs}). Default is NA.}
  \item{output.format}{character; Format of the output, \code{data.frame} (default), \link[sf]{sf-class}, \link[sftime]{sftime-class}, \link[terra]{SpatVector-class}, or \link[terra]{SpatRaster-class}.}
  \item{cpus}{numeric; Number of processing units. Default is detectCores()-1.}
  \item{progress}{logical; If progress bar is shown. Default is TRUE.}
  \item{soil3d}{logical; If 3D soil modellig is performed and \link{near.obs.soil} function is used for finding n nearest observations and distances to them. In this case, z position of the \code{data.staid.x.y.z} points to the depth column.}
  \item{depth.range}{numeric; Depth range for location mid depth in which to search for nearest observations (see function \link{near.obs.soil}). It's in the mid depth units. Default is 0.1.}
  \item{no.obs}{character; Possible values are \code{increase} (default) and \code{exactly}. If set to \code{increase}, in case if there is no \code{n.obs} observations in \code{depth.range} for a specific location, the \code{depth.range} is increased (multiplied by 2, 3, ...) until the number of observations are larger or equal to \code{n.obs}. If set to \code{exactly}, the function will raise an error when it come to the first location with no \code{n.obs} observations in specified \code{depth.range} (see function \link{near.obs.soil}).}
  \item{...}{Further arguments passed to \link[ranger]{predict.ranger} function, such as \code{type = "quantile"} and \code{quantiles = c(0.1,0.5,0.9)} for quantile regression, etc.}
}

\value{
A \code{data.frame}, \link[sf]{sf-class}, \link[sftime]{sftime-class}, \link[terra]{SpatVector-class}, or \link[terra]{SpatRaster-class} object (depends on \code{output.format} argument) with prediction - \code{pred} or \code{quantile..X.X} (quantile regression) columns.
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
library(ranger)
library(sp)
library(sf)
library(terra)
library(meteo)

# prepare data
demo(meuse, echo=FALSE)
meuse <- meuse[complete.cases(meuse@data),]
data = st_as_sf(meuse, coords = c("x", "y"), crs = 28992, agr = "constant")
# data = terra::vect(meuse)
# data.frame
# data <- as.data.frame(meuse)
# data$id = 1:nrow(data)
# data.staid.x.y.z <- c("id","x","y",NA)
fm.RFSI <- as.formula("zinc ~ dist + soil + ffreq")

# fit the RFSI model
rfsi_model <- rfsi(formula = fm.RFSI,
                   data = data, # meuse.df (use data.staid.x.y.z)
                   # data.staid.x.y.z = data.staid.x.y.z, # only if class(data) == data.frame
                   zero.tol = 0,
                   n.obs = 5, # number of nearest observations
                   # s.crs = st_crs(data), # nedded only if the coordinates are lon/lat (WGS84)
                   # p.crs = st_crs(data), # nedded only if the coordinates are lon/lat (WGS84)
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
                   quantreg = FALSE)
                   # quantreg = TRUE) # for quantile regression

rfsi_model
# OOB prediction error (MSE):       47758.14 
# R squared (OOB):                  0.6435869 
sort(rfsi_model$variable.importance)
sum("obs" == substr(rfsi_model$forest$independent.variable.names, 1, 3))

# Make RFSI prediction
# data.frame
# newdata <- as.data.frame(meuse.grid)
# newdata$id <- 1:nrow(newdata)
# newdata <- meuse.grid
newdata <- terra::rast(meuse.grid)
class(newdata)

# prediction
rfsi_prediction <- pred.rfsi(model = rfsi_model,
                             data = data, # meuse.df (use data.staid.x.y.z)
                             obs.col = "zinc",
                             # data.staid.x.y.z = data.staid.x.y.z, # only if class(data) == data.frame
                             newdata = newdata, # meuse.grid.df (use newdata.staid.x.y.z)
                             # newdata.staid.x.y.z = c("id", "x", "y", NA), # only if class(newdata) == data.frame
                             output.format = "SpatRaster", # "sf", # "SpatVector", 
                             zero.tol = 0,
                             # s.crs = st_crs(data), # meuse@proj4string, # NA # st_crs(data)
                             # newdata.s.crs = st_crs(data), # meuse@proj4string, # NA
                             # p.crs = st_crs(data), # meuse@proj4string, # NA
                             cpus = 1, # detectCores()-1,
                             progress = TRUE,
                             # type = "quantiles", # for quantile regression
                             # quantiles = c(0.1, 0.5, 0.9) # for quantile regression
)
class(rfsi_prediction)
names(rfsi_prediction)
# head(rfsi_prediction)

plot(rfsi_prediction)
# plot(rfsi_prediction['pred'])
# plot(rfsi_prediction['quantile..0.1'])
# plot(rfsi_prediction['quantile..0.5'])
# plot(rfsi_prediction['quantile..0.9'])
}

