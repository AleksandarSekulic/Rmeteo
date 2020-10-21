\name{near.obs}

\alias{near.obs}

\title{Finds n nearest observations from given locations.}

\description{The function finds n nearest observations from given locations and creates an object of \link[base]{data.frame} class, where first n columns are Euclidean distances to n nearest locations and next n columns are observations at n nearest stations, and rows are given locations. Further more it can calculate averages in circles with different radiuses, can find nearest observation in quadrants and calculate IDW predictions from nearest observations. It is based on \link[nabor]{knn} function of package nabor.}

\usage{
near.obs(locations,
         locations.x.y = c(1,2),
         observations,
         observations.x.y = c(1,2),
         zcol = 3,
         n.obs = 10,
         rm.dupl = TRUE,
         avg = FALSE,
         increment,
         range,
         direct = FALSE,
         idw=FALSE,
         idw.p=2)
}

\arguments{
  \item{locations}{data.frame with x and y coordinates columns, or \link[sp]{SpatialPoints}, \link[sp]{SpatialPointsDataFrame} or \link[sp]{SpatialPixelsDataFrame} object. Locations from which distances to observations are calculated.}
    \item{locations.x.y}{numeric or character vector; Positions or names of the x and y columns in \code{locations} if data.frame. Default is c(1,2).}
  \item{observations}{data.frame with x, y and observation columns, or \link[sp]{SpatialPoints} or \link[sp]{SpatialPointsDataFrame} object with observation column. Observations to which distances to locations are calculated.}
  \item{observations.x.y}{numeric or character vector; positions or names of the x and y columns in \code{observations} if data.frame. Default is c(1,2).}
  \item{zcol}{numeric or character; Column name or number showing the position of observations in \code{observations}. Default is 3.}
  \item{n.obs}{numeric; Number of nearest observations to be found. Note that it cannot be larger than number of obsevrations. Default is 10.}
  \item{rm.dupl}{boolean; Will duplicates, i.e. nearest observations where Euclidean distance is 0, be removed from the result. Default is TRUE.}
  \item{avg}{boolean; Will averages in circles with different radiuses be calculated. Default is FALSE.}
  \item{increment}{numeric; Increment of radiuses for calculation of averages in circles with different radiuses. Units depends on CRS of coordinates.}
  \item{range}{numeric; Maximum radius for calculation of averages in circles with different radiuses. Units depends on CRS of coordinates.}
  \item{direct}{boolean; Will nearest observation in quadrants be calculated. Default is FALSE.}
  \item{idw}{boolean; Will IDW predictions from \code{n.obs} nearest observations be calculated. Default is FALSE.}
  \item{idw.p}{numeric; Exponent parameter for IDW weights. Default is 2.}
}

\value{ \link[base]{data.frame} object.  Rows represents given locations. First \code{n.obs} columns are Euclidean distances to \code{n.obs} nearest observations. Next \code{n.obs} columns are observations at \code{n.obs} nearest stations. The following columns are averages in circles with different radiuses if \code{avg} is set to True. The following columns are nearest observation in quadrants if \code{direct} is set to True. The following columns are IDW prediction from nearest observation if \code{idw} is set to True.
}

\author{ALeksandar Sekulic \email{asekulic@grf.bg.ac.rs}}

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
# prepare data
# load observation - data.frame of mean temperatures
data(dtempc) 
str(dtempc)
data(stations)
#
str(stations)
lonmin=18 ;lonmax=22.5 ; latmin=40 ;latmax=46
library(sp)
library(spacetime)
serbia = point.in.polygon(stations$lon, stations$lat, c(lonmin,lonmax,lonmax,lonmin), 
                          c(latmin,latmin,latmax,latmax))
stations= stations[ serbia!=0, ] # stations in Serbia approx.
# create STFDF
temp <- meteo2STFDF(dtempc,stations, crs= CRS('+proj=longlat +datum=WGS84'))
str(temp)
# get one day (SpatialPointsDataFrame)
temp <- temp[, 1]
# find 5 nearest observations and distances to them (remove duplicates)
nearest_obs <- near.obs(locations = temp,
                        observations = temp,
                        zcol = "tempc",
                        n.obs = 5,
                        rm.dupl = TRUE)
str(nearest_obs)
}
