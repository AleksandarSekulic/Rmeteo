\name{pred.strk}

\alias{pred.strk}

\title{Spatio-temporal regression kriging prediction}

\description{Function for spatio-temporal regression kriging prediction based on \link[gstat]{krigeST}. The prediction is made for raster objects, i.e. for \link[spacetime]{STFDF-class}, \link[spacetime]{STSDF-class}, or for data.frame objects.}

\usage{
pred.strk(data,
          zcol=1,
          data.staid.x.y.time = c(1,2,3,4),
          obs,
          obs.staid.time = c(1,2),
          stations,
          stations.staid.x.y = c(1,2,3),
          newdata,
          newdata.staid.x.y.time = c(1,2,3),
          zero.tol=0,
          reg.coef,
          vgm.model,
          sp.nmax,
          time.nmax,
          by = 'time',
          tiling= FALSE,
          ntiles=64,
          output.format = "stfdf",
          parallel.processing = FALSE,
          pp.type = "snowfall",
          cpus=detectCores()-1,
          computeVar=FALSE,
          progress=TRUE,
          ...)
}

\arguments{
  \item{data}{\link[spacetime]{STFDF-class}, \link[spacetime]{STSDF-class} or \code{data.frame}; Contains dependent variable (observations) and covariates in space and time used to perform STRK. If \code{data.frame} object, it should have next columns: station ID (staid), longitude (x), latitude (y), time of the observation (time), observation value (obs) and covariates (cov1, cov2, ...). Covariate names should be the same as in the \code{reg.coef} (see below). If covariates are missing, overlay with \code{newdata} is tried. If overlay with \code{newdata} is not possible, then spatio-temporal ordinary kriging is performed.}
  \item{zcol}{numeric or character; Column name or number showing the position of dependent variable (observations) in \code{data}. Default is 1.}
  \item{data.staid.x.y.time}{numeric or character vector; Positions or names of the station ID (staid), longitude (x), latitude (y) and time columns in \code{data} if \code{data} is data.frame. Default is c(1,2,3,4).}
  \item{obs}{\code{data.frame}; Contains dependent variable (observations) and covariates in space and time. It should have next columns: station ID (staid), time of the observation (time), observation value (obs) and covariates (cov1, cov2, ...). This object is used together with \code{stations} (see below) to create \link[spacetime]{STFDF-class} object (if \code{data} object is missing) which is then used to perform STRK.  If covariates are missing, spatio-temporal ordinary kriging is performed.}
  \item{obs.staid.time}{numeric or character vector; Positions or names of the station ID (staid) and time columns in \code{obs}. Default is c(1,2).}
  \item{stations}{\code{data.frame}; It should have next columns: station ID (staid), longitude (x) and latitude (y) of the stations. This object is used together with \code{obs} (see above) if \code{data} object is missing.}
  \item{stations.staid.x.y}{numeric or character vector; Positions or names of the station ID (staid), longitude (x) and latitude (y) columns in \code{stations}. Default is c(1,2,3).}
  \item{newdata}{\link[spacetime]{STFDF-class}, \link[spacetime]{STSDF-class} or \code{data.frame}; Contains prediction locations and covariates. If \code{data.frame} object, it should have next columns: prediction location ID (staid), longitude (x), latitude (y), time of the prediction (time) and covariates (cov1, cov2, ...). Covariate names should be the same as in the \code{reg.coef} (see below).}
  \item{newdata.staid.x.y.time}{numeric or character vector; Positions or names of the prediction location ID (staid), longitude (x), latitude (y) and time columns in \code{newdata}. Default is c(1,2,3).}
  \item{zero.tol}{numeric; A distance value below (or equal to) which locations are considered as duplicates. Default is 0. See \link{rm.dupl}. Duplicates are removed to avoid singular covariance matrices in kriging.}
 \item{reg.coef}{named numeric vector; Named linear regression coefficients. Names of the coefficients (e.g. "Intercept", "temp_geo", "modis", "dem", "twi") will be used to match appropriate covariates from \code{data} (or \code{obs}) and \code{newdata}. Coefficients for metorological variables (temperature, precipitation, etc.) can be taken from data(\link{tregcoef}) or can be specified by the user.}
  \item{vgm.model}{StVariogramModel list; Spatio-temporal variogram of regression residuals (or observations if spatio-temporal ordinary kriging). See \link[gstat]{vgmST}. Spatio-temporal variogram model on residuals for metorological variables (temperature, precipitation, etc.) can be taken from data(\link{tvgms}) or can be specified by the user as a \link[gstat]{vgmST} object.}
  \item{sp.nmax}{numeric; A number of spatially nearest observations that should be used for kriging predictions. If /code{tiling} is TRUE (see below), then is a number of spatially nearest observations that should be used for each tile. Deafult is 20.}
  \item{time.nmax}{numeric; A number of temporally nearest observations that should be used for kriging predictions Deafult is 2.}
  \item{by}{cahracter; Will foreach loop by time (default) or station. If station is set, \code{sp.nmax} will be used for each station prediction.}
  \item{tiling}{logical; Should simplified local kriging be used. Default is FALSE. If TRUE, area is divided in tiles and kriging calculation is done for each tile separately. Number of observation used per tile is defined with \code{sp.nmax} and \code{time.nmax}. If \code{FALSE}, temporal local kriging will be applied defined with \code{time.nmax}, and \code{sp.nmax} will be ignored.}
  \item{ntiles}{numeric; A number of tiles for tilling. Default is 64. Ideally, each tile should contain less observations than \code{sp.nmax} and observations fall in neighboring tiles.}
  \item{output.format}{character; Format of the output, "STFDF" (default), "STSDF" or "data.frame" (\code{data.frame}).}
  \item{parallel.processing}{logical; If parallel processing is performed. Default is FALSE.}
  \item{pp.type}{character; Type (R package) of parallel processing, "snowfall" (default) or "doParallel".}
  \item{cpus}{numeric; Number of processing units. Default is detectCores()-1.}
  \item{computeVar}{logical; If kriging variance is computed. Ddefault is FALSE.}
  \item{progress}{logical; If progress bar is shown. Default is TRUE.}
  \item{...}{Further arguments passed to \link[gstat]{krigeST}.}
}

\value{
A \link[spacetime]{STFDF-class}, \link[spacetime]{STSDF-class} or \code{data.frame} object (depends on \code{output.format} argument), with columns:
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
library(gstat)
library(plyr)
# prepare data
# load observation - data.frame of mean temperatures
data(dtempc)
data(stations)

serbia= point.in.polygon(stations$lon, stations$lat, c(18,22.5,22.5,18), c(40,40,46,46))
st= stations[ serbia!=0, ]
dtempc <- dtempc[dtempc$staid \%in\% st$staid, ]
dtempc <- dtempc[complete.cases(dtempc),]

# create STFDF
stfdf <- meteo2STFDF(dtempc,st)
# Adding CRS
stfdf@sp@proj4string <- CRS('+proj=longlat +datum=WGS84')

# load covariates for mean temperatures
data(regdata)
# str(regdata)
regdata@sp@proj4string <- CRS('+proj=longlat +datum=WGS84')

# Calculate prediction of mean temperatures for "2011-07-06" 
# global model is used for regression and variogram
# load precalculated variograms
data(tvgms)
data(tregcoef)
results <- pred.strk(data = stfdf,
                     newdata = regdata[,2,drop=FALSE],
                     # newdata = regdata[,,drop=FALSE],
                     output.format = "stfdf",
                     reg.coef = tregcoef[[1]],
                     vgm.model = tvgms[[1]],
                     sp.nmax = 20,
                     time.nmax = 2,
                     computeVar=TRUE
)

# plot prediction
stplot(results[,,"pred", drop=F], col.regions=bpy.colors())
# stplot(results[,,"var", drop=F], col.regions=bpy.colors())

# Example with data.frames and parallel processing
library(snowfall)
library(doParallel)
# create data.frames
stfdf.df <- join(dtempc, st)
regdata.df <- as.data.frame(regdata)

results <- pred.strk(data = stfdf.df,
                     zcol=3,
                     data.staid.x.y.time = c(1,4,5,2),
                     # obs = stfdf.df, # if used, comment data argument
                     # obs.staid.time = c(1,2),
                     # stations = stfdf.df,
                     # stations.staid.x.y = c(1,4,5),
                     newdata = regdata.df[regdata.df$time=="2011-07-06", ],
                     # newdata = regdata.df,
                     newdata.staid.x.y.time = c(3,1,2,4),
                     reg.coef = tregcoef[[1]],
                     vgm.model = tvgms[[1]],
                     sp.nmax = 20,
                     time.nmax = 2,
                     parallel.processing = TRUE,
                     pp.type = "snowfall", # "doParallel"
                     cpus=detectCores()-1,
                     computeVar=TRUE,
                     progress=TRUE
)

# plot prediction
stplot(results[,,"pred", drop=F], col.regions=bpy.colors())
# stplot(results[,,"var", drop=F], col.regions=bpy.colors())

}

