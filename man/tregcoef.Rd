\name{tregcoef}
\docType{data}
\encoding{utf8}
\alias{tregcoef}
\title{Multiple linear regression coefficients for global and local daily air temperatures}
\description{
Multiple linear regression coefficients for mean, minimum, maximum daily temperature on geometric temperature trend (GTT), MODIS LST, digital elevation model (DEM) and topographic wetness index (TWI). The models are computed from GSOD, ECA&D, GHCN-Daily or local meteorological stations.
}
\usage{data(tregcoef)}
\format{
  A list of 9 multiple linear regression coefficients for daily air temperatures.
  \describe{
 \item{\code{tmeanGSODECAD}}{Multiple linear regression coefficients of global mean daily temperature on GTT, MODIS LST, DEM and TWI. Data used: GSOD and ECA&D}
\item{\code{tmeanGSODECAD_noMODIS}}{Multiple linear regression coefficients of global mean daily temperature on GTT, DEM, and TWI. Data used: GSOD and ECA&D}
\item{\code{tminGSODECAD}}{Multiple linear regression coefficients of global minimum daily temperature on GTT, MODIS LST, DEM, and TWI. Data used: GSOD and ECA&D}
\item{\code{tminGHCND}}{Multiple linear regression coefficients of global minimum daily temperature on GTT, MODIS LST, DEM, and TWI. Data used: GHCN-Daily}
\item{\code{tminGSODECAD_noMODIS}}{Multiple linear regression coefficients of global minimum daily temperature on GTT, DEM, and TWI. Data used: GSOD and ECA&D}
\item{\code{tmaxGSODECAD}}{Multiple linear regression coefficients of global maximum daily temperature on GTT, MODIS LST, DEM, and TWI. Data used: GSOD and ECA&D}
\item{\code{tmaxGHCND}}{Multiple linear regression coefficients of global maximum daily temperature on GTT, MODIS LST, DEM, and TWI. Data used: GHCN-Daily}
\item{\code{tmaxGSODECAD_noMODIS}}{Multiple linear regression coefficients of global maximum daily temperature on GTT, DEM, and TWI. Data used: GSOD and ECA&D}
\item{\code{tmeanHR}}{Multiple linear regression coefficients of Croatian mean daily temperature on GTT, DEM, and TWI. Data used: Croatian mean daily temperature dataset for the year 2008 (Croatian national meteorological network)}
   
  }
}
\author{Milan Kilibarda \email{kili@grf.bg.ac.rs}, Aleksandar Sekulic \email{asekulic@grf.bg.ac.rs}}
\references{
Kilibarda, M., T. Hengl, G. B. M. Heuvelink, B. Graeler, E. Pebesma, M. Percec Tadic, and B. Bajat (2014), Spatio-temporal interpolation of daily temperatures for global land areas at 1 km resolution, J. Geophys. Res. Atmos., 119, 2294-2313, doi:10.1002/2013JD020803.
}
\examples{
data(tregcoef)
tregcoef[[1]] # model for mean daily temp.
}
\keyword{datasets}
