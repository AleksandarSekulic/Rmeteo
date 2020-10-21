\name{tvgms}
\docType{data}
\encoding{utf8}
\alias{tvgms}
\title{Spatio-temporal variogram models for global and local daily air temperatures}
\description{
Variograms of residuals from multiple linear regression of mean, minimum, maximum daily temperatures on geometric temperature trend (GTT), MODIS LST, digital elevation model (DEM) and topographic wetness index (TWI). The models is computed from GSOD, ECA&D, GHCN-Daily or local meteorological stations. The obtained global or local models for mean, minimum, and maximum temperature can be used to produce gridded images of daily temperatures at high spatial and temporal resolution.
}
\usage{data(tvgms)}
\format{
  A list of 9 space-time sum-metric models for daily air temperatures, units: space km, time days.
  \describe{
 \item{\code{tmeanGSODECAD}}{Variogram for residuals from multiple linear regression of global mean daily temperature on GTT, MODIS LST, DEM, and TWI. D used: GSOD and ECA&D}
\item{\code{tmeanGSODECAD_noMODIS}}{Variogram for residuals from multiple linear regression of global mean daily temperature on GTT, DEM, and TWI. Data used: GSOD and ECA&D}
\item{\code{tminGSODECAD}}{Variogram for residuals from multiple linear regression of global minimum daily temperature on GTT, MODIS LST, DEM, and TWI. Data used: GSOD and ECA&D}
\item{\code{tminGHCND}}{Variogram for residuals from multiple linear regression of global minimum daily temperature on GTT, MODIS LST, DEM, and TWI. Data used: GHCN-Daily}
\item{\code{tminGSODECAD_noMODIS}}{Variogram for residuals from multiple linear regression of global minimum daily temperature on GTT, DEM, and TWI. Data used: GSOD and ECA&D}
\item{\code{tmaxGSODECAD}}{Variogram for residuals from multiple linear regression of global maximum daily temperature on GTT, MODIS LST, DEM, and TWI. Data used: GSOD and ECA&D}
\item{\code{tmaxGHCND}}{Variogram for residuals from multiple linear regression of global maximum daily temperature on GTT, MODIS LST, DEM, and TWI. Data used: GHCN-Daily}
\item{\code{tmaxGSODECAD_noMODIS}}{Variogram for residuals from multiple linear regression of global maximum daily temperature on GTT,  DEM, and TWI. Data used: GSOD and ECA&D}
\item{\code{tmeanHR}}{Variogram for residuals from multiple linear regression of Croatian mean daily temperature on GTT, DEM, and TWI. Data used: Croatian mean daily temperature dataset for the year 2008 (Croatian national meteorological network)}
   
  }
}
\author{Milan Kilibarda \email{kili@grf.bg.ac.rs}, ALeksandar Sekulic \email{asekulic@grf.bg.ac.rs}}
\references{
Kilibarda, M., T. Hengl, G. B. M. Heuvelink, B. Graeler, E. Pebesma, M. Percec Tadic, and B. Bajat (2014), Spatio-temporal interpolation of daily temperatures for global land areas at 1 km resolution, J. Geophys. Res. Atmos., 119, 2294-2313, doi:10.1002/2013JD020803.
}
\examples{
data(tvgms)
tvgms[[1]] # model for mean daily temp.
}
\keyword{datasets}
