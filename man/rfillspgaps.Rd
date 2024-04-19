\name{rfillspgaps}
\alias{rfillspgaps}
\title{Close gaps of a grid or raster Layer data}
\description{The function close gaps of a raster data by using IDW.}
\usage{
rfillspgaps(rasterLayer,
            maskPol=NULL,
            nmax =50,
            zcol=1,
            ...)
}
\arguments{
\item{rasterLayer}{\link[terra]{SpatRaster-class}, \link[raster]{RasterLayer-class}, \link[sp]{SpatialGrid-class} or \link[sp]{SpatialPixels-class}; Raster that contains NAs.}
\item{maskPol}{\link[sf:sf]{sf-class}, \link[terra]{SpatVector-class}, \link[sp]{SpatialPolygons} or \link[sp]{SpatialPolygonsDataFrame}; Area of interest to spatially fill \code{rasterLayer} missing values and to mask \code{rasterLayer}.}
\item{nmax}{see \link[gstat]{krige}, idw function.}
\item{zcol}{\code{integer} or \code{character}; variable column name or number showing position of a variable in \code{rasterLayer} to be interpolated.}
\item{...}{arguments passed to \link[gstat]{krige}, idw function.}
}
\value{raster object with NA replaced using IDW in \code{rasterLayer} format.}
\author{Milan Kilibarda \email{kili@grf.bg.ac.rs}}
\references{
Kilibarda, M., T. Hengl, G. B. M. Heuvelink, B. Graeler, E. Pebesma, M. Percec Tadic, and B. Bajat (2014), Spatio-temporal interpolation of daily temperatures for global land areas at 1 km resolution, J. Geophys. Res. Atmos., 119, 2294-2313, doi:10.1002/2013JD020803;

Kilibarda M., M. Percec Tadic, T. Hengl, J. Lukovic, B. Bajat - Spatial Statistics (2015), Global geographic and feature space coverage of temperature data in the context of spatio-temporal interpolation, doi:10.1016/j.spasta.2015.04.005.
}
\seealso{
\code{\link{rfilltimegaps}}
\code{\link{pred.strk}}
}
\examples{
   library(terra)
   data(nlmodis20110712)
   data(NLpol)
   
   nlmodis20110712 <- terra::rast(nlmodis20110712)
   # SpaVector
   NLpol = vect(NLpol)
   crs(NLpol) <- "epsg:4326"
   # # sf
   # NLpol <- st_as_sf(NLpol) #, crs = st_crs(4326))
   
   \donttest{
   plot(nlmodis20110712)
   }
   # fill spatial gaps
   r=rfillspgaps(nlmodis20110712,NLpol)
   \donttest{
   plot(r)
   }
  
}
