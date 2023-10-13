\name{dtempc_ogimet}
\docType{data}
\encoding{utf8}
\alias{dtempc_ogimet}
\title{Mean daily temperature in degrees Celsius for the year 2019 for Serbia}
\description{Sample data set of mean daily temperature measurements from the OGIMET service for the year 2019 for Serbian territory.}
\usage{data(dtempc_ogimet)}
\format{
The \code{dtempc_ogimet} contains the following columns:
  \describe{
  \item{\code{staid}}{character; station ID from OGIMET}
  \item{\code{name}}{character; station name}
  \item{\code{lon}}{numeric; Longitude}
  \item{\code{lat}}{numeric; Latitude}
  \item{\code{elevation}}{numeric; Hight}
  \item{\code{time}}{Date; day of the measurement}
  \item{\code{tmean}}{numeric; mean daily temperature in degrees Celsius}
  \item{\code{dem}}{numeric; Digital Elevation Model (DEM) in meters}
  \item{\code{twi}}{numeric; Topographic Wetness Index (TWI)}
  \item{\code{cdate}}{numeric; Cumulative day from 1960}
  \item{\code{doy}}{numeric; Day of year}
  \item{\code{gtt}}{numeric; Geometrical temperature trend}
}
}
\author{Aleksandar Sekulic \email{asekulic@grf.bg.ac.rs}}
\references{
\itemize{
\item OGIMET service (\url{https://www.ogimet.com/home.phtml.en})
}
}
\examples{
## load data 
data(dtempc_ogimet)
str(dtempc)
}
\keyword{datasets}