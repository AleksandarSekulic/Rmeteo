\name{pred.strk}

\alias{pred.strk}

\title{Spatio-temporal regression kriging prediction}

\description{Function for spatio-temporal regression kriging prediction based on \link[gstat]{krigeST}.}

\usage{
pred.strk(data,
          obs.col=1,
          data.staid.x.y.z = NULL,
          newdata,
          newdata.staid.x.y.z = NULL,
          z.value = NULL,
          crs = NA,
          zero.tol=0,
          reg.coef,
          vgm.model,
          sp.nmax=20,
          time.nmax=2,
          by='time',
          tiling= FALSE,
          ntiles=64,
          output.format = "STFDF",
          parallel.processing = FALSE,
          pp.type = "snowfall",
          cpus=detectCores()-1,
          computeVar=FALSE,
          progress=TRUE,
          ...)
}

\arguments{
  \item{data}{\link[spacetime]{STFDF-class}, \link[spacetime]{STSDF-class}, \link[spacetime]{STIDF-class}, \link[sf:sf]{sf-class}, \link[sftime:st_sftime]{sftime-class}, \link[terra]{SpatVector-class} or \link[base]{data.frame}; Contains target variable (observations) and covariates in space and time used to perform STRK. If \link[base]{data.frame} object, it should have next columns: station ID (staid), longitude (x), latitude (y), 3rd component - time, depth, ... (z) of the observation, observation value (obs), and covariates (cov1, cov2, ...). Covariate names should be the same as in the \code{reg.coef} (see below). If covariates are missing, overlay with \code{newdata} is tried. If overlay with \code{newdata} is not possible, then spatio-temporal ordinary kriging is performed.}
  \item{obs.col}{numeric or character; Column name or number showing position of the observation column in the \code{data}. Default is 1.}
  \item{data.staid.x.y.z}{numeric or character vector; Positions or names of the station ID (staid), longitude (x), latitude (y) and 3rd component - time, depth (z) columns in \link[base]{data.frame} object (e.g. c(1,2,3,4)). If \code{data} is \link[sf:sf]{sf-class}, \link[sftime:st_sftime]{sftime-class}, or \link[terra]{SpatVector-class} object, \code{data.staid.x.y.z} is used to point staid and z position. Default is NULL.}
  \item{newdata}{\link[spacetime]{STFDF-class}, \link[spacetime]{STSDF-class}, \link[spacetime]{STIDF-class}, \link[sf:sf]{sf-class}, \link[sftime:st_sftime]{sftime-class}, \link[terra]{SpatVector-class}, \link[terra]{SpatRaster-class} or \link[base]{data.frame}; Contains prediction locations and covariates used for STRK prediction. If \link[base]{data.frame} object, it should have next columns: prediction location ID (staid), longitude (x), latitude (y), 3rd component - time, depth, ... (z), and covariates (cov1, cov2, ...). Covariate names have to be the same as in the \code{reg.coef} (see below).}
  \item{newdata.staid.x.y.z}{numeric or character vector; Positions or names of the prediction location ID (staid), longitude (x), latitude (y) and 3rd component (z) columns in data.frame \code{newdata} object (e.g. c(1,2,3,4)). If \code{newdata} is \link[sf:sf]{sf-class}, \link[sftime:st_sftime]{sftime-class}, \link[terra]{SpatVector-class} or \link[terra]{SpatRaster-class} object, \code{newdata.staid.x.y.z} is used to point staid and z position. Default is NULL.}
  \item{z.value}{vector; A vector of 3rd component - time, depth, ... (z) values if \code{newdata} is \link[terra]{SpatRaster-class}.}
  \item{crs}{\link[sf]{st_crs} or \link[terra]{crs}; Source CRS of \code{data} and \code{newdata}. If \code{data} or \code{newdata} contains crs, \code{crs} will not be used. Default is NA.}
  \item{zero.tol}{numeric; A distance value below (or equal to) which locations are considered as duplicates. Default is 0. See \link{rm.dupl}. Duplicates are removed to avoid singular covariance matrices in kriging.}
   \item{reg.coef}{numeric; Vector of named linear regression coefficients. Names of the coefficients (e.g. "Intercept", "temp_geo", "modis", "dem", "twi") will be used to match appropriate covariates from \code{data}. Coefficients for metorological variables (temperature, precipitation, etc.) can be taken from data(\link{tregcoef}) or can be specified by the user.}
  \item{vgm.model}{StVariogramModel list; Spatio-temporal variogram of regression residuals (or observations if spatio-temporal ordinary kriging). See \link[gstat]{vgmST}. Spatio-temporal variogram model on residuals for metorological variables (temperature, precipitation, etc.) can be taken from data(\link{tvgms}) or can be specified by the user as a \link[gstat]{vgmST} object.}
  \item{sp.nmax}{numeric; A number of spatially nearest observations that should be used for kriging predictions. If \code{tiling} is TRUE (see below), then is a number of spatially nearest observations that should be used for each tile. Deafult is 20.}
  \item{time.nmax}{numeric; A number of temporally nearest observations that should be used for kriging predictions Deafult is 2.}
  \item{by}{cahracter; Will foreach loop by time (default) or station. If station is set, \code{sp.nmax} will be used for each station prediction.}
  \item{tiling}{logical; Should simplified local kriging be used. Default is FALSE. If TRUE, area is divided in tiles and kriging calculation is done for each tile separately. Number of observation used per tile is defined with \code{sp.nmax} and \code{time.nmax}. If \code{FALSE}, temporal local kriging will be applied defined.}
  \item{ntiles}{numeric; A number of tiles for tilling. Default is 64. Ideally, each tile should contain less observations than \code{sp.nmax} and observations fall in neighboring tiles.}
  \item{output.format}{character; Format of the output, \link[spacetime]{STFDF-class} (default), \link[spacetime]{STSDF-class}, \link[spacetime]{STIDF-class}, \link[base]{data.frame}, \link[sf:sf]{sf-class}, \link[sftime:st_sftime]{sftime-class}, \link[terra]{SpatVector-class}, or \link[terra]{SpatRaster-class}.}
  \item{parallel.processing}{logical; If parallel processing is performed. Default is FALSE.}
  \item{pp.type}{character; Type (\code{R} package) for parallel processing, "snowfall" (default) or "doParallel".}
  \item{cpus}{numeric; Number of processing units. Default is detectCores()-1.}
  \item{computeVar}{logical; If kriging variance is computed. Default is FALSE.}
  \item{progress}{logical; If progress bar is shown. Default is TRUE.}
  \item{...}{Further arguments passed to \link[gstat]{krigeST} function.}
}

\value{
A \link[spacetime]{STFDF-class}, \link[spacetime]{STSDF-class}, \link[spacetime]{STIDF-class}, \link[base]{data.frame}, \link[sf:sf]{sf-class}, \link[sftime:st_sftime]{sftime-class}, \link[terra]{SpatVector-class}, or \link[terra]{SpatRaster-class} object (depends on \code{output.format} argument), with columns (elements):
  \item{pred}{Predictions.}
  \item{tlm}{Trend.}
  \item{var}{Kriging variance, if \code{computeVar=TRUE}.}
}

\author{Milan Kilibarda \email{kili@grf.bg.ac.rs}, Aleksandar Sekulic \email{asekulic@grf.bg.ac.rs}}

\references{
Kilibarda, M., T. Hengl, G. B. M. Heuvelink, B. Graeler, E. Pebesma, M. Percec Tadic, and B. Bajat (2014), Spatio-temporal interpolation of daily temperatures for global land areas at 1 km resolution, J. Geophys. Res. Atmos., 119, 2294-2313, doi:10.1002/2013JD020803.
}

\seealso{
 \code{\link{tregcoef}}
 \code{\link{tvgms}}
 \code{\link{regdata}}
 \code{\link{meteo2STFDF}}
 \code{\link{tgeom2STFDF}}
}

\examples{
library(sp)
library(spacetime)
library(sf)
library(gstat)
library(plyr)
# prepare data
# load observation - data.frame of mean temperatures
# preparing data
data(dtempc) 
data(stations)
data(regdata) # covariates, made by mete2STFDF function
regdata@sp@proj4string <- CRS('+proj=longlat +datum=WGS84')

lonmin=18 ;lonmax=22.5 ; latmin=40 ;latmax=46
serbia = point.in.polygon(stations$lon, stations$lat, c(lonmin,lonmax,lonmax,lonmin), 
                          c(latmin,latmin,latmax,latmax))
st = stations[ serbia!=0, ] # stations in Serbia approx.
obs.staid.time = c("staid", "time")
stations.staid.lon.lat = c(1,2,3)
crs = CRS('+proj=longlat +datum=WGS84')
delta = NULL

\dontrun{
# create STFDF
stfdf <- meteo2STFDF(obs = dtempc,
                     stations = st,
                     crs = crs)

# Calculate prediction of mean temperatures for "2011-07-05" and "2011-07-06" 
# global model is used for regression and variogram
# load precalculated variograms
data(tvgms) # ST variogram models
data(tregcoef) # MLR coefficients

### Example with STFDF and without parallel processing
results <- pred.strk(data = stfdf, # observations
                     newdata = regdata, # prediction locations with covariates
                     # newdata = regdata[,2,drop=FALSE], # for one day only
                     output.format = "STFDF", # data.frame | sf | sftime | SpatVector | SpatRaster 
                     reg.coef = tregcoef[[1]], # MLR coefficients
                     vgm.model = tvgms[[1]], # STRK variogram model
                     sp.nmax = 20,
                     time.nmax = 2,
                     computeVar=TRUE
)
class(results)
# plot prediction
# results@sp=as(results@sp,'SpatialPixelsDataFrame')
# stplot(results[,,"pred", drop= FALSE], col.regions=bpy.colors())
# stplot(results[,,"var", drop= FALSE], col.regions=bpy.colors())
}

\dontrun{
### Example with data.frames and parallel processing - SpatRaster output
library(terra)
library(doParallel)
# create data.frame
stfdf.df <- join(dtempc, st)
summary(stfdf.df)
regdata.df <- as.data.frame(regdata)
results <- pred.strk(data = stfdf.df,
                     obs.col = 3,
                     data.staid.x.y.z = c(1,4,5,2),
                     newdata = regdata.df,
                     newdata.staid.x.y.z = c(3,1,2,4),
                     crs = CRS("EPSG:4326"),
                     output.format = "SpatRaster", # STFDF |data.frame | sf | sftime | SpatVector
                     reg.coef = tregcoef[[1]],
                     vgm.model = tvgms[[1]],
                     sp.nmax = 20,
                     time.nmax = 2,
                     parallel.processing = TRUE,
                     pp.type = "doParallel", # "snowfall"
                     cpus = 1, # detectCores()-1,
                     computeVar = TRUE,
                     progress = TRUE
)
# plot prediction
# plot(results$`2011-07-06`[["pred"]])
# plot(results$`2011-07-06`[["var"]])
}
}

