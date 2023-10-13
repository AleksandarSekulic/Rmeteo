\name{data.prepare}

\alias{data.prepare}

\title{Prepare data}

\description{Function for data preparation for RFSI and STRK functions. It transforms data to a data.frame.}

\usage{
data.prepare(data,
             data.staid.x.y.z=NULL,
             obs.col=NULL,
             s.crs=NA
)
}

\arguments{
  \item{data}{\link[sf:sf]{sf-class}, \link[sftime:st_sftime]{sftime-class}, \link[terra]{SpatVector-class}, \link[terra]{SpatRaster-class} or \link[base]{data.frame}; Contains target variable (observations) and covariates. If \link[base]{data.frame} object, it should have next columns: station ID (staid), longitude (x), latitude (y), 3rd component - time, depth, ... (z) of the observation, and observation value (obs).}
  \item{data.staid.x.y.z}{numeric or character vector; Positions or names of the station ID (staid), longitude (x), latitude (y) and 3rd component (z) columns in \link[base]{data.frame} object (e.g. c(1,2,3,4)). If \code{data} is \link[sf:sf]{sf-class}, \link[sftime:st_sftime]{sftime-class}, or \link[terra]{SpatVector-class} object, \code{data.staid.x.y.z} is used to point staid and z position. Set z position to NA (e.g. c(1,2,3,NA)) or ommit it (e.g. c(1,2,3)) for spatial interpolation. Default is NULL.}
  \item{obs.col}{numeric or character; Column name or number showing position of the observation column in the \code{data}. Default is 1.}
  \item{s.crs}{\link[sf]{st_crs} or \link[terra]{crs}; Source CRS of \code{data}. If \code{data} contains crs, \code{s.crs} will not be used. Default is NA.}
}

\value{
A list with the following elements:
  \item{data.df}{A \link[base]{data.frame} obtained from \code{data}.}
  \item{data.staid.x.y.z}{Positions of the station ID (staid), longitude (x), latitude (y) and 3rd component (z) columns in \link[base]{data.frame} object (e.g. c(1,2,3,4)).}
  \item{s.crs}{Source CRS of \code{data}.}
  \item{obs.col}{Column number showing position of the observation column in the \code{data}.}
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
library(sf)
library(meteo)
library(sp)

# prepare data
demo(meuse, echo=FALSE)
meuse <- meuse[complete.cases(meuse@data),]
data = st_as_sf(meuse, coords = c("x", "y"), crs = 28992, agr = "constant")
data.df <- data.prepare(data,
                        obs.col="zinc")
str(data.df)
}

