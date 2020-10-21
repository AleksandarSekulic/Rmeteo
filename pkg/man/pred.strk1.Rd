\name{pred.strk1}
\alias{pred.strk1}
\title{Spatio-temporal regression kriging}
\description{ Function for spatio-temporal regression kriging prediction based on \link[gstat]{krigeST}. The prediction is made for points given in \link[base]{data.frame}.}
\usage{
pred.strk1(obs, stations, newdata, zero.tol=0,
reg.coef=list(
tmean=c(-0.126504415,0.4051734447,0.4943247727,0.0001837527,-0.0189207588),
tmin = c(-0.9825601517,0.5672140021,0.3344561638, 0.0003119777,-0.0243629638),
tmax = c(1.7873573081,0.350228076, 0.5569091092, 0.0002571338,-0.0012988123)
)[['tmean']],
vgm.model=list(tmean=vgmST("sumMetric",
                                    space=vgm( 14.13, "Sph", 5903, 1.933),
                                    time =vgm(0, "Sph",  0.1, 0),
                                    joint=vgm(9.06, "Sph", 2054, 0.474),
                                    stAni=497.9),
                         tmin = vgmST("sumMetric",
                                    space=vgm( 22.682, "Sph", 5725, 3.695),
                                    time =vgm(0, "Sph",  0.1, 0),
                                    joint=vgm(9.457, "Sph",1888, 1.67),
                                    stAni=485),
                         tmax = vgmST("sumMetric",
                                    space=vgm( 8.31, "Sph", 4930, 2.872),
                                    time =vgm(0, "Sph",  0.1, 0),
                                    joint=vgm(11.175, "Sph", 2117, 1.75),
                                    stAni=527) ) [['tmean']] ,
computeVar=FALSE, out.remove=FALSE, threshold.res=15, returnList=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs}{ data.frame; observations data frame contains station ID column, time column (day of observation), measured variable column and ovariates columns (in exactly that order). It can contain additional variables (columns).}
  \item{stations}{ data.frame; Stations data frame contains at least station ID column, longitude (or x) and latitude (or y) column (in exactly that order).It can contain additional variables (columns).}
  \item{newdata}{ data.frame; predictions data frame contains longitude (or x) and latitude (or y) column, time column (day of prediction) and covariates columns (in exactly that order)}
  \item{zero.tol}{distance values less than or equal to this threshold value locations are considered as duplicates, see \link{rm.dupl}, duplicates are removed to avoid singular covariance matrices in kriging.}
  \item{reg.coef}{linear regression coefficients; order is assumed as intercept, \code{dynamic.cov}, \code{static.cov}. Coefficients can be specified by user; depending on type, number and order of dynamic and static covariates.  At the moment the function contains regression coefficient for mean, minimum and maximum temperature calculated globally for GSOD and ECA&D data set on geometrical temperature trend, MODIS LST-8 day, elevation and TWI, see \link{regdata}. Coefficients for mean temperature are defined by default.}
  \item{vgm.model}{spatio-temporal variogram of regression residuals, see \link[gstat]{vgmST}. At the moment the function contains spatio-temporal variogram model on residuals for mean, minimum and maximum temperature calculated globally for GSOD and ECA&D data set. Regression residuals on geometrical temperature trend, MODIS LST-8 day, elevation and TWI, see \link{regdata}. Ranges are in km. Spatio-temporal variogram for mean temperatures is defined by default. User can specified own variogram model as \link[gstat]{vgmST} object.}
  \item{computeVar}{if TRUE, just variance is computed}
  \item{out.remove}{if TRUE, potential outliers are removed. Removing procedure is iterative, all location with residual higher than defined threshold (\code{treshold.res}) are selected. Only location with highest cross validation residual is removed, than cross validation is done again, the procedure removing one by one location run until all locations have residuals under defined threshold.  
}
  \item{threshold.res}{critical threshold for removing potential outliers}
  \item{returnList}{if TRUE, result islist; if FALSE, result is data frame }
}

\value{
An list object containing:
  \item{pred}{an object of \link[sp]{SpatialPointsDataFrame-class} with column contains prediction or variance}
  \item{cv}{cross validation information for points used in prediction, as object of \link[spacetime]{STFDF-class}}
  \item{remst}{removed locations as an object of \link[sp]{Spatial-class}, if \code{out.remove=TRUE} }
  \item{remobs}{removed locations with observations as an object of \link[spacetime]{STFDF-class}, if \code{out.remove=TRUE} }
A data frame object contains longitude (or x) and latitude (or y) column, time column (day of prediction) and prediction value.
}

\author{Aleksandar Sekulic \email{asekulic@grf.bg.ac.rs}, Milan Kilibarda \email{kili@grf.bg.ac.rs}  }
\references{
Kilibarda, M., T. Hengl, G. B. M. Heuvelink, B. Graeler, E. Pebesma, M. Percec Tadic, and B. Bajat (2014), Spatio-temporal interpolation of daily temperatures for global land areas at 1 km resolution, J. Geophys. Res. Atmos., 119, 2294-2313, doi:10.1002/2013JD020803.
}
\seealso{
 \code{\link{regdata}}
 \code{\link{meteo2STFDF}}
 \code{\link{tgeom2STFDF}}
 
}
\examples{
## prepare data
## load observation - data.frame of mean temperatures
data(dtempc) 
str(dtempc)
data(stations)
library(sp)
library(spacetime)
library(gstat)

#str(stations)
## lonmin,lonmax,lonmax, lonmin   latmin, latnmin,latmax,latmax
serbia= point.in.polygon(stations$lon, stations$lat, c(18,22.5,22.5,18), c(40,40,46,46))
st= stations[ serbia!=0, ]
## create STFDF
temp <- meteo2STFDF(dtempc,st)
rm(dtempc)
# str(temp)
## Adding CRS
temp@sp@proj4string <- CRS('+proj=longlat +datum=WGS84')

## load covariates for mean temperatures
data(regdata)
# str(regdata)
regdata@sp@proj4string <- CRS('+proj=longlat +datum=WGS84')

#creating newdata
pred <- data.frame(regdata@sp@coords[1:5,1],regdata@sp@coords[1:5,2], 
                  '2011-07-05', regdata@data$temp_geo[1:5], regdata@data$modis[1:5], 
                  regdata@sp@data$dem[1:5], regdata@sp@data$twi[1:5])
## pred names
names(pred)=c("x", "y", "time", "temp_geo", "modis", "dem", "twi")

#creating observation
gg <- regdata
time <- gg@time
gg@data$dem = rep(gg@sp@data[,1],length(time))
gg@data$twi = rep(gg@sp@data[,2],length(time))

temp_geo <- sapply(1:length(time),
                function(i) over(temp@sp,as(gg[,i,'temp_geo'],'SpatialPixelsDataFrame')))
modis <- sapply(1:length(time),
                function(i) over(temp@sp,as(gg[,i,'modis'],'SpatialPixelsDataFrame' ) ) )
dem <- sapply(1:length(time),
              function(i) over(temp@sp,as(gg[,i,'dem'],'SpatialPixelsDataFrame' ) ) )
twi <- sapply(1:length(time),
                function(i) over(temp@sp,as(gg[,i,'twi'],'SpatialPixelsDataFrame' ) ) )

temp_geo <- do.call('cbind',temp_geo)
temp_geo <- as.vector(temp_geo)
modis <- do.call('cbind',modis)
modis <- as.vector(modis)
dem <- do.call('cbind',dem)
dem <- as.vector(dem)
twi <- do.call('cbind',twi)
twi <- as.vector(twi)

t1 <- which(as.character(index(time[1])) == as.character( index(temp@time)) )
t2 <- which(as.character(index(time[ length(time) ])) == as.character(index(temp@time)) )

temp <- temp[,t1:t2, drop=FALSE]

temp$temp_geo = temp_geo
temp$modis = modis
temp$dem = dem
temp$twi = twi

obs = as.data.frame(temp)[,c(7,4,11,12,13,14,15)]

# Calculate prediction of points in pred data frame
# global model is used for regression and variogram
# load precalculated variograms
data(tvgms)
data(tregcoef)
res= pred.strk1(obs,st, newdata= pred, 
              reg.coef=tregcoef[[1]] ,vgm.model=tvgms[[1]], returnList=TRUE )

#str(res)

# res1= pred.strk1(obs,st, newdata= pred, reg.coef=tregcoef[[1]] ,vgm.model=tvgms[[1]], 
#               returnList=FALSE, out.remove=TRUE, threshold.res=15 )

# str(res1)

}
