\name{stations_ogimet}
\docType{data}
\encoding{utf8}
\alias{stations_ogimet}
\title{Data frame containing stations' information from the OGIMET service for Serbian territory}
\description{Data frame containing stations' information of daily observations from the OGIMET service for the year 2019 for Serbian territory.}
\usage{data(stations_ogimet)}
\format{
The \code{dtempc_ogimet} contains the following columns:
  \describe{
  \item{\code{staid}}{character; station ID from OGIMET}
  \item{\code{name}}{character; station name}
  \item{\code{lon}}{numeric; Longitude}
  \item{\code{lat}}{numeric; Latitude}
  \item{\code{elevation}}{numeric; Hight}
  \item{\code{dem}}{numeric; Digital Elevation Model (DEM) in meters}
  \item{\code{twi}}{numeric; Topographic Wetness Index (TWI)}
}
}

\author{Aleksandar Sekulic \email{asekulic@grf.bg.ac.rs}}
\references{
\itemize{
\item OGIMET service (\url{https://www.ogimet.com/})
}
}
\examples{
# load data:
data(stations_ogimet)
str(stations)
library(sp)
coordinates(stations) <-~ lon +lat
stations@proj4string <-CRS('+proj=longlat +datum=WGS84')
\donttest{
plot(stations)
}
}
\keyword{datasets}
