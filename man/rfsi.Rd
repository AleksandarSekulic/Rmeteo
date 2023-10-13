\name{rfsi}

\alias{rfsi}

\title{Random Forest Spatial Interpolation (RFSI) model}

\description{Function for creation of Random Forest Spatial Interpolation (RFSI) model (Sekulić et al. 2020). Besides environmental covariates, RFSI uses additional spatial covariates: (1) observations at n nearest locations and (2) distances to them, in order to include spatial context into the random forest.}

\usage{
rfsi(formula,
     data,
     data.staid.x.y.z = NULL,
     n.obs = 5,
     avg = FALSE,
     increment = 10000,
     range = 50000,
     quadrant = FALSE,
     use.idw = FALSE,
     idw.p = 2,
     s.crs = NA,
     p.crs = NA,
     cpus = detectCores()-1,
     progress = TRUE,
     soil3d = FALSE,
     depth.range = 0.1,
     no.obs = 'increase',
     ...)
}

\arguments{
  \item{formula}{formula; Formula for specifying target variable and covariates (without nearest observations and distances to them). If \code{z~1}, an RFSI model using only nearest obsevrations and distances to them as covariates will be made.}
  \item{data}{\link[sf:sf]{sf-class}, \link[sftime:st_sftime]{sftime-class}, \link[terra]{SpatVector-class} or \link[base]{data.frame}; Contains target variable (observations) and covariates used for making an RFSI model. If \link[base]{data.frame} object, it should have next columns: station ID (staid), longitude (x), latitude (y), 3rd component - time, depth, ... (z) of the observation, observation value (obs) and covariates (cov1, cov2, ...). If covariates are missing, the RFSI model using only nearest obsevrations and distances to them as covariates (\code{formula=z~1}) will be made.}
  \item{data.staid.x.y.z}{numeric or character vector; Positions or names of the station ID (staid), longitude (x), latitude (y) and 3rd component (z) columns in \link[base]{data.frame} object (e.g. c(1,2,3,4)). If \code{data} is \link[sf:sf]{sf-class}, \link[sftime:st_sftime]{sftime-class}, or \link[terra]{SpatVector-class} object, \code{data.staid.x.y.z} is used to point staid and z position. Set z position to NA (e.g. c(1,2,3,NA)) or ommit it (e.g. c(1,2,3)) for spatial interpolation. Default is NULL.}
  \item{n.obs}{numeric; Number of nearest observations to be used as covariates in RFSI model (see function \link{near.obs}). Note that it cannot be larger than number of obsevrations. Default is 5.}
  \item{avg}{boolean; Averages in circles covariate - will averages in circles with different radiuses be calculated (see function \link{near.obs}). Default is FALSE.}
  \item{increment}{numeric; Increment of radiuses for calculation of averages in circles with different radiuses (see function \link{near.obs}). Units depends on CRS.}
  \item{range}{numeric; Maximum radius for calculation of averages in circles with different radiuses (see function \link{near.obs}). Units depends on CRS.}
  \item{quadrant}{boolean; Nearest observations in quadrants covariate - will nearest observation in quadrants be calculated (see function \link{near.obs}). Default is FALSE.}
  \item{use.idw}{boolean; IDW prediction as covariate - will IDW predictions from \code{n.obs} nearest observations be calculated (see function \link{near.obs}). Default is FALSE.}
  \item{idw.p}{numeric; Exponent parameter for IDW weights (see function \link{near.obs}). Default is 2.}
  \item{s.crs}{\link[sf]{st_crs} or \link[terra]{crs}; Source CRS of \code{data}. If \code{data} contains crs, \code{s.crs} will be overwritten. Default is NA.}
  \item{p.crs}{\link[sf]{st_crs} or \link[terra]{crs}; Projection CRS for \code{data} reprojection. If NA, \code{s.crs} will be used for distance calculation. Note that observations should be in projection for finding nearest observations based on Eucleadean distances (see function \link{near.obs}). Default is NA.}
  \item{cpus}{numeric; Number of processing units. Default is detectCores()-1.}
  \item{progress}{logical; If progress bar is shown. Default is TRUE.}
  \item{soil3d}{logical; If 3D soil modellig is performed and \link{near.obs.soil} function is used for finding n nearest observations and distances to them. In this case, z position of the \code{data.staid.x.y.z} points to the depth column.}
  \item{depth.range}{numeric; Depth range for location mid depth in which to search for nearest observations (see function \link{near.obs.soil}). It's in the mid depth units. Default is 0.1.}
  \item{no.obs}{character; Possible values are \code{increase} (default) and \code{exactly}. If set to \code{increase}, in case if there is no \code{n.obs} observations in \code{depth.range} for a specific location, the \code{depth.range} is increased (multiplied by 2, 3, ...) until the number of observations are larger or equal to \code{n.obs}. If set to \code{exactly}, the function will raise an error when it come to the first location with no \code{n.obs} observations in specified \code{depth.range} (see function \link{near.obs.soil}).}
  \item{...}{Further arguments passed to \link[ranger]{ranger}, such as \code{quantreg}, \code{importance}, etc.}
}

\value{
RFSI model of class \link[ranger]{ranger}.
}

\note{Observations should be in projection for finding nearest observations based on Eucleadean distances (see function \link{near.obs}). If crs is not specified in the \code{data} object or through the \code{s.crs} parameter, the coordinates will be used as they are in projection. Use \code{s.crs} and \code{p.crs} if the coordinates of the \code{data} object are in lon/lat (WGS84).}

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
# data.staid.x.y.z <- c(15,"x","y",NA)
fm.RFSI <- as.formula("zinc ~ dist + soil + ffreq")

# fit the RFSI model
rfsi_model <- rfsi(formula = fm.RFSI,
                   data = data, # meuse.df (use data.staid.x.y.z)
                   # data.staid.x.y.z = data.staid.x.y.z, # only if class(data) == data.frame
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

rfsi_model
# OOB prediction error (MSE):       47758.14 
# R squared (OOB):                  0.6435869 
sort(rfsi_model$variable.importance)
sum("obs" == substr(rfsi_model$forest$independent.variable.names, 1, 3))

}

