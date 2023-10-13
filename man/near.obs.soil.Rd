\name{near.obs.soil}

\alias{near.obs.soil}

\title{Finds n nearest observations from given locations for soil mapping.}

\description{The function finds n nearest observations from given locations and at specific depth range from location min depth and creates an object of \link[base]{data.frame} class.
First n columns are Euclidean distances to n nearest locations and next n columns are observations at n nearest stations, and rows are given locations. It is based on \link[nabor]{knn} function of package nabor.}

\usage{
near.obs.soil(locations,
              locations.x.y.md = c(1,2,3),
              observations,
              observations.x.y.md = c(1,2,3),
              obs.col = 4,
              n.obs = 5,
              depth.range = 0.1,
              no.obs = 'increase',
              parallel.processing = TRUE,
              pp.type = "doParallel", # "snowfall"
              cpus = detectCores()-1)
}

\arguments{
  \item{locations}{data.frame with x and y coordinates and mid depth columns, or \link[sf:sf]{sf-class}, \link[terra]{SpatVector-class} or \link[terra]{SpatRaster-class} object; Locations (FROM) for which n nearest observations are found and distances are calculated.}
  \item{locations.x.y.md}{numeric or character vector; Positions or names of the x, y, and mid depth columns in \code{locations} if data.frame. Default is c(1,2,3).}
  \item{observations}{data.frame with x, y, mid depth and observation columns, or \link[sf:sf]{sf-class} or \link[terra]{SpatVector-class} object with mid depth and observation columns; Observations (TO).}
  \item{observations.x.y.md}{numeric or character vector; positions or names of the x, y, and mid depth columns in \code{observations} if data.frame. Default is c(1,2,3).}
  \item{obs.col}{numeric or character; Column name or number showing position of the observation column in the \code{observations}. Default is 4.}
  \item{n.obs}{numeric; Number of nearest observations to be found. Note that it cannot be larger than number of obsevrations. Default is 5.}
  \item{depth.range}{numeric; Depth range for location mid depth in which to search for nearest observations. It's in the mid depth units. Default is 0.1.}
  \item{no.obs}{character; Possible values are \code{increase} (default) and \code{exactly}. If set to \code{increase}, in case if there is no \code{n.obs} observations in \code{depth.range} for a specific location, the \code{depth.range} is increased (multiplied by 2, 3, ...) until the number of observations are larger or equal to \code{n.obs}. If set to \code{exactly}, the function will raise an error when it come to the first location with no \code{n.obs} observations in specified \code{depth.range}.}
  \item{parallel.processing}{logical; If parallel processing is performed. Default is FALSE.}
  \item{pp.type}{character; Type (R package) used for parallel processing, "doParallel" (default) or "snowfall".}
  \item{cpus}{numeric; Number of processing units. Default is detectCores()-1.}
}

\value{\link[base]{data.frame} object.  Rows represents specific locations. First \code{n.obs} columns are Euclidean distances to \code{n.obs} nearest observations. Next \code{n.obs} columns are observations at \code{n.obs} nearest stations.
}

\author{Aleksandar Sekulic \email{asekulic@grf.bg.ac.rs}, Anatol Helfenstein \email{anatol.helfenstein@wur.nl}}

\note{The function is intended for soil mapping applications.}

\seealso{
  \code{\link[nabor]{knn}}
  \code{\link{near.obs}}
  \code{\link{rfsi}}
  \code{\link{pred.rfsi}}
  \code{\link{tune.rfsi}}
  \code{\link{cv.rfsi}}
}

\examples{
library(sp)
library(sf)
library(meteo)
# prepare data
# load observation - data.frame of mean temperatures
demo(meuse, echo=FALSE)
meuse <- meuse[complete.cases(meuse@data),]
locations = st_as_sf(meuse, coords = c("x", "y"), crs = 28992, agr = "constant")
locations = # terra::rast(meuse.grid)
observations = st_as_sf(meuse, coords = c("x", "y"), crs = 28992, agr = "constant")
# find 5 nearest observations and distances to them (remove duplicates)
nearest_obs <- near.obs.soil(locations = locations, # from
                             locations.x.y.md =  c("x","y","dist"),
                             observations = observations, # to
                             observations.x.y.md= c("x","y","dist"),
                             obs.col = "zinc",
                             n.obs = 5) # number of nearest observations
str(nearest_obs)
summary(nearest_obs)
}
