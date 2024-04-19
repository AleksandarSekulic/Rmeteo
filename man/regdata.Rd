\name{regdata}
\alias{regdata}
\docType{data}
\encoding{utf8}
\alias{regdata}
\title{Dynamic and static covariates for spatio-temporal regression kriging}
\description{Dynamic and static covariates for spatio-temporal regression kriging of \link[spacetime]{STFDF-class}. The \code{regdata} contains geometrical temperature trend, MODIS LST 8-day splined at daily resolution, elevation and topographic wetness index. }
\usage{data(regdata)}
\format{
The \code{regdata} contains the following dynamic and static covariates:
  \describe{
  \item{\code{regdata$temp_geo}}{numeric; geometrical temperature trend for mean temperature, calculated with \code{\link{tgeom2STFDF}} ; from 2011-07-05 to 2011-07-09, in degree Celsius  }
  \item{\code{regdata$modis}}{numeric;  MODIS LST 8-day splined at daily resolution, missing pixels are filtered by spatial splines and 8-day values are splined at daily level; from 2011-07-05 to 2011-07-09, in degree Celsius  }
  \item{\code{regdata@sp$dem}}{numeric; elevation data obtained from Worldgrids (depricated)}
   \item{\code{regdata@sp$twi}}{numeric; SAGA Topographic Wetness Index (TWI) from Worldgrids (depricated)}
}
}
\author{Milan Kilibarda \email{kili@grf.bg.ac.rs}}
\examples{
data(regdata)
str(regdata)
library(sp) # spplot
library(spacetime) # stplot

\donttest{
stplot(regdata[,,'modis']) # plot modis data
spplot(regdata@sp,zcol='twi', col.regions = bpy.colors() ) # plot TWI
spplot(regdata@sp,zcol='dem', col.regions = bpy.colors() ) # plot dem
}
}
\keyword{datasets}
