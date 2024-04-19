\name{dem_twi_srb}
\docType{data}
\encoding{utf8}
\alias{dem_twi_srb}
\title{Digital Elevation Model (DEM) and Topographic Wetness Index (TWI) for Serbia}
\description{Digital Elevation Model (DEM) and Topographic Wetness Index (TWI) for Serbia in \link[terra]{SpatRaster} format.}
\usage{data(dem_twi_srb)}
\format{
The \code{dem_twi_srb} contains the following layers:
  \describe{
  \item{\code{dem}}{Digital Elevation Model (DEM) in meters}
  \item{\code{twi}}{Topographic Wetness Index (TWI)}
}
}
\author{Aleksandar Sekulic \email{asekulic@grf.bg.ac.rs}}
\examples{
library(terra)
# load data 
data(dem_twi_srb)
terra::unwrap(dem_twi_srb)
}
\keyword{datasets}
