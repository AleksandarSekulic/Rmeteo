\name{rfilltimegaps}
\alias{rfilltimegaps}
\title{Disaggregation in the time dimension through the use of splines for each pixel}
\description{The function creates an object of \link[spacetime]{STFDF-class} class, spatio-temporal data with full space-time grid, from another \link[spacetime]{STFDF-class} and fills attribute data for missing values in time using splines.}
\usage{
rfilltimegaps(stfdf,
              tunits="day",
              attrname=1,
              ...)
}
\arguments{
\item{stfdf}{\link[spacetime]{STFDF-class}; object with time information of minimum length 2, and gap in time dimension.}
\item{tunits}{\code{character}; increment of the sequence used to generete time infromation for temporal gap. See 'Details'.}
\item{attrname}{\code{integer} or \code{character}; varible from \link[spacetime]{STFDF-class} to be splined in time.}
\item{...}{arguments passed to \link[stats]{splinefun}, function spline.}
}
\details{
\code{tunits} can be specified in several ways:
\itemize{
\item{A number, taken to be in seconds}
\item{A object of class \link[base]{difftime}}
\item{A character string, containing one of \code{"sec"}, \code{"min"}, \code{"hour"}, \code{"day"}, \code{"DSTday"}, \code{"week"}, \code{"month"}, \code{"quarter"} or \code{"year"}. This can optionally be preceded by a (positive or negative) integer and a space, or followed by \code{"s"}}
}
The difference between \code{"day"} and \code{"DSTday"} is that the former ignores changes to/from daylight savings time and the latter takes the same clock time each day. (\code{"week"} ignores DST (it is a period of 144 hours), but \code{"7 DSTdays"}) can be used as an alternative. \code{"month"} and \code{"year"} allow for DST.)
}
\value{\link[spacetime]{STFDF-class} object with filled temporal gaps.}
\author{Milan Kilibarda \email{kili@grf.bg.ac.rs}}
\references{
Kilibarda, M., T. Hengl, G. B. M. Heuvelink, B. Graeler, E. Pebesma, M. Percec Tadic, and B. Bajat (2014), Spatio-temporal interpolation of daily temperatures for global land areas at 1 km resolution, J. Geophys. Res. Atmos., 119, 2294-2313, doi:10.1002/2013JD020803;

Kilibarda M., M. Percec Tadic, T. Hengl, J. Lukovic, B. Bajat - Spatial Statistics (2015), Global geographic and feature space coverage of temperature data in the context of spatio-temporal interpolation, doi:10.1016/j.spasta.2015.04.005.
}
\seealso{
\code{\link{rfillspgaps}}
\code{\link{pred.strk}}
}
\examples{
  data(nlmodis20110704)
  data(nlmodis20110712)
  data(NLpol)

  # fill spatial gaps
  library(raster)
  NLpol@proj4string <- nlmodis20110704@proj4string
  
  \donttest{
  nlmodis20110704 <- rfillspgaps(nlmodis20110704,NLpol)
  nlmodis20110712 <- rfillspgaps(nlmodis20110712,NLpol)
  
  nlmodis20110704 <- as(nlmodis20110704,"SpatialPixelsDataFrame")
  names(nlmodis20110704)='m1'
  nlmodis20110712 <- as(nlmodis20110712,"SpatialPixelsDataFrame")
  names(nlmodis20110712)='m2'
  
  nlmodis20110704@data <- cbind(nlmodis20110704@data, nlmodis20110712@data)
  
  df<-reshape(nlmodis20110704@data , varying=list(1:2), v.names="modis",direction="long", 
            times=as.Date(c('2011-07-04','2011-07-12')), ids=1:dim(nlmodis20110704)[1])
  
  library(spacetime)
  stMODIS<- STFDF(as( nlmodis20110704, "SpatialPixels"), 
                  time= as.Date(c('2011-07-04','2011-07-12')), 
                  data.frame(modis=df[,'modis']))
  
  stplot(stMODIS, col.regions=bpy.colors())
  stMODIS <- rfilltimegaps(stMODIS)
  stplot(stMODIS, col.regions=bpy.colors())
  }
}

