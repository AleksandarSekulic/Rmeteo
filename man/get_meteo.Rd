\name{get_meteo}
\alias{get_meteo}
\title{Get daily meteorological data for specific location and dates.}
\description{The function gives back daily metorological data for specific location and dates for Europe for 1991-2020 period (from MeteoEurope1km dataset).}
\usage{
get_meteo(loc,
          dates,
          var = "tmean",
          source = "MeteoEurope1km")
}
\arguments{
\item{loc}{\code{sf}, \code{SpatVector}, \link[base]{data.frame}, \code{matrix}, \code{numeric} or \code{integer}; Locations in WGS84 (EPSG:4326). If \link[base]{data.frame} or \code{matrix} columns are lon/lat.}
\item{dates}{\code{Date} or \code{character}; Date(s).}
\item{var}{\code{character}; Daily meteorological variable. Possible values are: \code{"tmean"} (default), \code{"tmax"}, \code{"tmin"}, \code{"prcp"}, or \code{"slp"}.}
\item{source}{\code{character}; Data source. Possible values are: \code{"MeteoEurope1km"} (default).}
}

\value{\link[base]{data.frame} object with daily meteorological values for specific locations (rows) and specific dates (columns). First column is location ID.}
\author{Aleksandar Sekulić \email{asekulic@grf.bg.ac.rs}}

\seealso{
\code{\link{pred.strk}}
}
\examples{
library(terra)

loc <- c(21, 45)
# loc <- as.data.frame(rbind(c(21, 45),
#                            c(21,45.5),
#                            c(21.5,45),
#                            c(21.5,45.5)))

dates <- as.Date("2020-12-25")
# dates <- seq(as.Date("2020-12-25"), as.Date("2020-12-31"), by="day")

tmean <- get_meteo(loc,
                   dates,
                   var = "tmean", # "tmax" "tmin" "prcp" "slp"
                   source = "MeteoEurope1km")
}

