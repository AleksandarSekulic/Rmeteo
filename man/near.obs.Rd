\name{near.obs}

\alias{near.obs}

\title{Finds n nearest observations from given locations.}

\description{The function finds n nearest observations from given locations and creates an object of \link[base]{data.frame} class.
First n columns are Euclidean distances to n nearest locations and next n columns are observations at n nearest stations, and rows are given locations. Further more it can calculate averages in circles with different radiuses, can find nearest observation in quadrants (directions) and calculate IDW predictions from nearest observations. It is based on \link[nabor]{knn} function of package nabor.}

\usage{
near.obs(locations,
         locations.x.y = c(1,2),
         observations,
         observations.x.y = c(1,2),
         obs.col = 3,
         n.obs = 10,
         rm.dupl = TRUE,
         avg = FALSE,
         increment,
         range,
         quadrant = FALSE,
         idw=FALSE,
         idw.p=2)
}

\arguments{
  \item{locations}{data.frame with x and y coordinates columns, or \link[sf:sf]{sf-class}, \link[terra]{SpatVector-class} or \link[terra]{SpatRaster-class} object; Locations (FROM) for which n nearest observations are found and distances are calculated.}
  \item{locations.x.y}{numeric or character vector; Positions or names of the x and y columns in \code{locations} if data.frame. Default is c(1,2).}
  \item{observations}{data.frame with x, y and observation columns, or \link[sf:sf]{sf-class} or \link[terra]{SpatVector-class} object with an observation column; Observations (TO).}
  \item{observations.x.y}{numeric or character vector; Positions or names of the x and y columns in \code{observations} if data.frame. Default is c(1,2).}
  \item{obs.col}{numeric or character; Column name or number showing position of the observation column in the \code{observations}. Default is 3.}
  \item{n.obs}{numeric; Number of nearest observations to be found. Note that it cannot be larger than number of obsevrations. Default is 10.}
  \item{rm.dupl}{boolean; Remove spatial duplicates - will the spatial duplicates (nearest observations where Euclidean distance is 0) be removed from the result. Default is TRUE.}
  \item{avg}{boolean; Averages in circles - will averages in circles with different radiuses be calculated. Default is FALSE.}
  \item{increment}{numeric; Increment of radiuses for calculation of averages in circles with different radiuses. Units depends on CRS.}
  \item{range}{numeric; Maximum radius for calculation of averages in circles with different radiuses. Units depends on CRS.}
  \item{quadrant}{boolean; Nearest observations in quadrants - will nearest observation in quadrants be calculated. Default is FALSE.}
  \item{idw}{boolean; IDW prediction as predictor - will IDW predictions from \code{n.obs} nearest observations be calculated. Default is FALSE.}
  \item{idw.p}{numeric; Exponent parameter for IDW weights. Default is 2.}
}

\value{ \link[base]{data.frame} object.  Rows represents specific locations. First \code{n.obs} columns are Euclidean distances to \code{n.obs} nearest observations. Next \code{n.obs} columns are observations at \code{n.obs} nearest stations. The following columns are averages in circles with different radiuses if \code{avg} is set to TRUE. The following columns are nearest observation in quadrants if \code{direct} is set to TRUE. The following columns are IDW prediction from nearest observation if \code{idw} is set to TRUE.
}

\author{Aleksandar Sekulic \email{asekulic@grf.bg.ac.rs}}

\references{
Sekulić, A., Kilibarda, M., Heuvelink, G. B., Nikolić, M. & Bajat, B. Random Forest Spatial Interpolation. Remote. Sens. 12, 1687, https://doi.org/10.3390/rs12101687 (2020).
}

\note{The function can be used in any case if it is needed to find n nearest observations from given locations and distances to them.}

\seealso{
  \code{\link[nabor]{knn}}
  \code{\link{rfsi}}
  \code{\link{pred.rfsi}}
  \code{\link{tune.rfsi}}
  \code{\link{cv.rfsi}}
}

\examples{
library(sp)
library(sf)
library(terra)
library(meteo)
# prepare data
# load observation - data.frame of mean temperatures
demo(meuse, echo=FALSE)
meuse <- meuse[complete.cases(meuse@data),]
locations = terra::rast(meuse.grid)
observations = st_as_sf(meuse, coords = c("x", "y"), crs = 28992, agr = "constant")
# find 5 nearest observations and distances to them (remove duplicates)
nearest_obs <- near.obs(locations = locations, # from
                        observations = observations, # to
                        obs.col = "zinc",
                        n.obs = 5, # number of nearest observations
                        rm.dupl = TRUE) 
str(nearest_obs)
summary(nearest_obs)
}
