\name{temp_geom}
\alias{temp_geom}
\title{Calculate geometrical temperature trend
}
\description{Calculate geometrical temperature trend for mean, minimum or maximum temperature.
}
\usage{
temp_geom(day,
          fi,
          variable="mean",
          ab = NULL)
}
\arguments{
  \item{day}{integer; Day of the year (from 1 to 366). Single value or vector of days of the year (only if \code{fi} is single value).}
  \item{fi}{numeric; Latitude. Single value or vector of latitudes (only if \code{day} is single value).}
  \item{variable}{character; Geometrical temperature trend calculated for mean, minimum or maximum temperature; Possible values are \code{'mean'}, \code{'min'} or \code{'max'}. \code{'mean'} is default.}
  \item{ab}{Predefined coefficients to be used instead of incorporated.}
}
\value{ 
A numerical single value or \code{vector} with calculated geometrical temperature trend.
}
\references{
Kilibarda, M., T. Hengl, G. B. M. Heuvelink, B. Graeler, E. Pebesma, M. Percec Tadic, and B. Bajat (2014), Spatio-temporal interpolation of daily temperatures for global land areas at 1 km resolution, J. Geophys. Res. Atmos., 119, 2294-2313, doi:10.1002/2013JD020803.
}
\author{Milan Kilibarda \email{kili@grf.bg.ac.rs}, Aleksandar Sekulic \email{asekulic@grf.bg.ac.rs}}
\examples{
tgeom <- temp_geom(day = 1,
          fi = 45.33,
          variable="mean")
          
tgeom_vect <- temp_geom(day = 1:365,
          fi = 45.33,
          variable="mean")
          
tgeom_vect2 <- temp_geom(day = 1,
          fi = seq(35, 45, 0.5),
          variable="mean")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

