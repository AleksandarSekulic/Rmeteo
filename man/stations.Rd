\name{stations}
\docType{data}
\encoding{utf8}
\alias{stations}
\title{Data frame containing stations' information }
\description{Data frame containing stations' information of merged daily observations from the Global Surface Summary of Day (GSOD) with European Climate Assessment & Data set (ECA&D) for the month July 2011.}
\usage{data(stations)}
\format{
The \code{stations} contains the following columns:
  \describe{
  \item{\code{staid}}{character; station ID from GSOD or ECA&D data set}
  \item{\code{lon}}{numeric; longitude coordinate}
  \item{\code{lat}}{numeric; longitude coordinate}
  \item{\code{elev_1m}}{numeric; elevation derived from station metadata in m}
  \item{\code{data_source}}{Factor; data source, GSOD or ECA&D}
  \item{\code{station_name}}{character; station name}
}
}
\author{Milan Kilibarda and Tomislav Hengl}
\references{
\itemize{
\item Global Surface Summary of the day data (\url{ftp://ftp.ncdc.noaa.gov/pub/data/gsod/}) 
\item European Climate Assessment & Dataset (\url{http://eca.knmi.nl/dailydata/predefinedseries.php}) 
}
}
\examples{
## load data:
data(stations)
str(stations)
library(sp)
coordinates(stations) <-~ lon +lat
stations@proj4string <-CRS('+proj=longlat +datum=WGS84')
plot(stations)
}
\keyword{datasets}
