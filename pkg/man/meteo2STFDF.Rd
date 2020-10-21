\name{meteo2STFDF}
\alias{meteo2STFDF}
\title{Create an object of \link[spacetime]{STFDF} class from two data frames (observation and stations) }
\description{The function creates an object of \link[spacetime]{STFDF} class, spatio-temporal data with full space-time grid, from two data frames (observation and stations). Observations data frame minimum contains station ID column, time column (day of observation) and measured variable column. Stations data frame contains at least station ID column, longitude (or x) and latitude (or y) column.}
\usage{meteo2STFDF(obs, stations, obs.staid.time = c(1, 2), 
            stations.staid.lon.lat = c(1, 2, 3), crs=CRS(as.character(NA)), delta=NULL )}
\arguments{
\item{obs}{ data.frame; observations data frame minimum contains station ID column, time column (day of observation) and measured variable column. It can contain additional variables (columns).}
  \item{stations}{ data.frame; Stations data frame contains at least station ID column, longitude (or x) and latitude (or y) column.It can contain additional variables (columns).}
  \item{obs.staid.time}{ numeric; records the column positions where in \code{obs} (observation) data frame the station ID and time values are stored.}
  \item{stations.staid.lon.lat}{ numeric; records the column positions where in \code{stations} data frame the station ID, longitude (x) and latitude (y) values are stored.}
  \item{crs}{ CRS; coordinate reference system (see \link[sp]{CRS}) of \code{stations} coordinates
}
  \item{delta}{ time;  time interval to end points in seconds
}

}
\value{ \link[spacetime]{STFDF} object}
\author{Milan Kilibarda \email{kili@grf.bg.ac.rs}, ALeksandar Sekulic \email{asekulic@grf.bg.ac.rs}}
\note{The function is intended for conversion of meteorological data to \link[spacetime]{STFDF} object, but can be used for similar spatio temporal data stored in two separated tables.}
\seealso{
\code{\link{tgeom2STFDF}}
\code{\link{pred.strk}}
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
st= stations[ serbia!=0, ] # stations in Serbia approx.
# create STFDF
temp <- meteo2STFDF(dtempc,st, crs= CRS('+proj=longlat +datum=WGS84'))
str(temp)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

