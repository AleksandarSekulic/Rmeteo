\name{near.obs.soil}

\alias{near.obs.soil}

\title{Nearest observatios for soil mapping.}

\description{The function finds n nearest observations from given locations and at specific depth range from location min depth and creates an object of \link[base]{data.frame} class, where first n columns are Euclidean distances to n nearest locations and next n columns are observations at n nearest stations, and rows are given locations. It is based on \link[nabor]{knn} function of package nabor.}

\usage{
near.obs.soil(locations,
              locations.x.y.md = c(1,2,3),
              observations,
              observations.x.y.md = c(1,2,3),
              zcol = 4,
              n.obs = 5,
              depth.range = 0.1,
              no.obs = 'increase')
}

\arguments{
  \item{locations}{data.frame with x and y coordinates and mid depth columns, or \link[sp]{SpatialPoints-class}, \link[sp]{SpatialPointsDataFrame-class} or \link[sp]{SpatialPixelsDataFrame-class} object. Locations from which distances to observations are calculated.}
  \item{locations.x.y.md}{numeric or character vector; Positions or names of the x, y, and mid depth columns in \code{locations} if data.frame. Default is c(1,2,3).}
  \item{observations}{data.frame with x, y, mid depth and observation columns, or \link[sp]{SpatialPoints-class} or \link[sp]{SpatialPointsDataFrame-class} object with mid depth and observation columns. Observations to which distances to locations are calculated.}
  \item{observations.x.y.md}{numeric or character vector; positions or names of the x, y, and mid depth columns in \code{observations} if data.frame. Default is c(1,2,3).}
  \item{zcol}{numeric or character; Column name or number showing the position of observations in \code{observations}. Default is 3.}
  \item{n.obs}{numeric; Number of nearest observations to be found. Note that it cannot be larger than number of observations. Default is 5.}
  \item{depth.range}{numeric; Depth range from location mid depth in which to search for nearest observations. It is in the mid depth units. Default is 0.1.}
  \item{no.obs}{character; Possible values are \code{increase} (default) and \code{exactly}. If set to \code{increase}, for locations where there is no \code{n.obs} observations in \code{depth.range}, the \code{depth.range} is increased (multiplied by 2, 3, ...) until there the number of observations are larger or equal to \code{n.obs}. If set to \code{exactly}, the function will raise an error for the first locations where there is no \code{n.obs} observations in \code{depth.range}.}
}

\value{ \link[base]{data.frame} object. Rows represents given locations. First \code{n.obs} columns are Euclidean distances to \code{n.obs} nearest observations. Next \code{n.obs} columns are observations at \code{n.obs} nearest stations.
}

\author{ALeksandar Sekulic \email{asekulic@grf.bg.ac.rs}, Anatol Helfenstein {anatol.helfenstein@wur.nl}}

\references{
}

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

}
