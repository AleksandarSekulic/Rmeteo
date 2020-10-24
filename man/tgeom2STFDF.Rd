\name{tgeom2STFDF}
\alias{tgeom2STFDF}
\title{Calculate geometrical temperature trend
}
\description{Calculate geometrical temperature trend for mean, minimum or maximum temperature.
}
\usage{
tgeom2STFDF(grid, time, variable = "mean",ab=NULL)
}
\arguments{
  \item{grid}{ object of \link[sp]{Spatial-class} (Points, Grid or Pixels) with associated coordinate reference systems (\link[sp]{CRS-class}). If                       CRS is not defined longitude latitude is assumed. }
  \item{time}{object holding time information, reasonably it is day (calendar date), or vector of days
}
  \item{variable}{ character;\code{'mean'}, \code{'min'} or \code{'max'} ; geometrical temperature trend is calculated for mean, minimum or maximum; \code{'mean'} is default. }
    \item{ab}{Predefined coefficients to be used instead of incorporated.}
}
\value{ \link[spacetime]{STFDF} object with calculated \code{temp_geo} geometrical temperature trend. The calculated values are stored in \code{obj@data} slot.
}
\references{
Kilibarda, M., T. Hengl, G. B. M. Heuvelink, B. Graeler, E. Pebesma, M. Percec Tadic, and B. Bajat (2014), Spatio-temporal interpolation of daily temperatures for global land areas at 1 km resolution, J. Geophys. Res. Atmos., 119, 2294-2313, doi:10.1002/2013JD020803.
}
\author{Milan Kilibarda \email{kili@grf.bg.ac.rs}, ALeksandar Sekulic \email{asekulic@grf.bg.ac.rs}}
\examples{
library(sp)
library(spacetime)
## create one point from lon lat
pos <- SpatialPoints(coords = cbind(19.22,45.33)) 
## temp_geom for 1st Jan 2011
tg1 <- tgeom2STFDF(pos,as.POSIXct("2011-01-01") )  
tg1

## temp_geom for the 2011 at pos location
tg365<- tgeom2STFDF(pos,time = seq(as.POSIXct("2011-01-01"), as.POSIXct("2011-12-31"), 
                    by="day") ) 
stplot(tg365, mode='ts')

data(regdata) 
## DEM and TWI data for Serbia at 1 km resolution
# str(regdata@sp)
spplot(regdata@sp, zcol='dem', col.regions=bpy.colors() )

## temp_geom for Serbia 1st and 2nd Jully 2011
tgSrb<- tgeom2STFDF(regdata@sp,time = seq(as.POSIXct("2011-07-01"), 
                    as.POSIXct("2011-07-02"), by="day") ) 

## temp_geom for "2011-07-01" , "2011-07-02"
# stplot(tgSrb, col.regions = bpy.colors() ) 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

