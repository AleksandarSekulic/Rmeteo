\name{tune.rfsi}

\alias{tune.rfsi}

\title{Tuning of Random Forest Spatial Interpolation (RFSI) model}

\description{Function for tuning of Random Forest Spatial Interpolation (RFSI) model using k-fold leave-location-out cross-validation (Sekulić et al. 2020).}

\usage{
tune.rfsi(formula,
          data,
          data.staid.x.y.z = NULL,
          use.idw = FALSE,
          s.crs = NA,
          p.crs = NA,
          tgrid,
          tgrid.n=10,
          tune.type = "LLO",
          k = 5,
          seed=42,
          folds,
          acc.metric,
          fit.final.model=TRUE,
          cpus = detectCores()-1,
          progress = TRUE,
          soil3d = FALSE,
          no.obs = 'increase',
          ...)
}

\arguments{
  \item{formula}{formula; Formula for specifying target variable and covariates (without nearest observations and distances to them). If \code{z~1}, an RFSI model using only nearest obsevrations and distances to them as covariates will be tuned.}
  \item{data}{\link[sf:sf]{sf-class}, \link[sftime:st_sftime]{sftime-class}, \link[terra]{SpatVector-class} or \link[base]{data.frame}; Contains target variable (observations) and covariates used for making an RFSI model. If \link[base]{data.frame} object, it should have next columns: station ID (staid), longitude (x), latitude (y), 3rd component - time, depth, ... (z) of the observation, observation value (obs) and covariates (cov1, cov2, ...). If covariates are missing, the RFSI model using only nearest obsevrations and distances to them as covariates (\code{formula=z~1}) will be tuned.}
  \item{data.staid.x.y.z}{numeric or character vector; Positions or names of the station ID (staid), longitude (x), latitude (y) and 3rd component (z) columns in \link[base]{data.frame} object (e.g. c(1,2,3,4)). If \code{data} is \link[sf:sf]{sf-class}, \link[sftime:st_sftime]{sftime-class}, or \link[terra]{SpatVector-class} object, \code{data.staid.x.y.z} is used to point staid and z position. Set z position to NA (e.g. c(1,2,3,NA)) or ommit it (e.g. c(1,2,3)) for spatial interpolation. Default is NULL.}
  \item{use.idw}{boolean; IDW prediction as covariate - will IDW predictions from \code{n.obs} nearest observations be calculated and tuned (see function \link{near.obs}). Default is FALSE.}
  \item{s.crs}{\link[sf]{st_crs} or \link[terra]{crs}; Source CRS of \code{data}. If \code{data} contains crs, \code{s.crs} will be overwritten. Default is NA.}
  \item{p.crs}{\link[sf]{st_crs} or \link[terra]{crs}; Projection CRS for \code{data} reprojection. If NA, \code{s.crs} will be used for distance calculation. Note that observations should be in projection for finding nearest observations based on Eucleadean distances (see function \link{near.obs}). Default is NA.}
  \item{tgrid}{data.frame; Possible tuning parameters. The column names are same as the tuning parameters. Possible tuning parameters are: \code{n.obs}, \code{num.trees}, \code{mtry}, \code{min.node.size}, \code{sample.fraction}, \code{splirule}, \code{idw.p}, and \code{depth.range}.}
  \item{tgrid.n}{numeric; Number of randomly chosen \code{tgrid} combinations for tuning of RFSI. If larger than \code{tgrid}, will be set to \code{length(tgrid)}}
  \item{tune.type}{character; Type of cross-validation: leave-location-out ("LLO"), leave-time-out ("LTO") - TO DO, and leave-location-time-out ("LLTO") - TO DO. Default is "LLO".}
  \item{k}{numeric; Number of random folds that will be created with \link[CAST]{CreateSpacetimeFolds} function if \code{folds} is column. Default is 5.}
  \item{seed}{numeric; Random seed that will be used to generate folds with \link[CAST]{CreateSpacetimeFolds} function.}
  \item{folds}{numeric or character vector or value; Showing folds column (if value) or rows (vector) of \code{data} observations used for cross-validation. If missing, will be created with \link[CAST]{CreateSpacetimeFolds} function.}
  \item{acc.metric}{character; Accuracy metric that will be used as a criteria for choosing an optimal RFSI model. Possible values for regression: "ME", "MAE", "NMAE", "RMSE" (default), "NRMSE", "R2", "CCC". Possible values for classification: "Accuracy","Kappa" (default), "AccuracyLower", "AccuracyUpper", "AccuracyNull", "AccuracyPValue", "McnemarPValue".}
  \item{fit.final.model}{boolean; Fit the final RFSI model. Defailt is TRUE.}
  \item{cpus}{numeric; Number of processing units. Default is detectCores()-1.}
  \item{progress}{logical; If progress bar is shown. Default is TRUE.}
  \item{soil3d}{logical; If 3D soil modellig is performed and \link{near.obs.soil} function is used for finding n nearest observations and distances to them. In this case, z position of the \code{data.staid.x.y.z} points to the depth column.}
  \item{no.obs}{character; Possible values are \code{increase} (default) and \code{exactly}. If set to \code{increase}, in case if there is no \code{n.obs} observations in \code{depth.range} for a specific location, the \code{depth.range} is increased (multiplied by 2, 3, ...) until the number of observations are larger or equal to \code{n.obs}. If set to \code{exactly}, the function will raise an error when it come to the first location with no \code{n.obs} observations in specified \code{depth.range} (see function \link{near.obs.soil}).}
  \item{...}{Further arguments passed to \link[ranger]{ranger}.}
}

\value{
A list with elements:
  \item{combinations}{data.frame; All tuned parameter combinations with chosen accuracy metric value.}
  \item{tuned.parameters}{numeric vector; Tuned parameters with chosen accuracy metric value.}
  \item{final.model}{\link[ranger]{ranger}; Final RFSI model (if \code{fit.final.model=TRUE}).}
}

\author{Aleksandar Sekulic \email{asekulic@grf.bg.ac.rs}}

\references{
Sekulić, A., Kilibarda, M., Heuvelink, G. B., Nikolić, M. & Bajat, B. Random Forest Spatial Interpolation. Remote. Sens. 12, 1687, https://doi.org/10.3390/rs12101687 (2020).
}

\seealso{
  \code{\link{near.obs}}
  \code{\link{rfsi}}
  \code{\link{pred.rfsi}}
  \code{\link{cv.rfsi}}
}

\examples{
library(CAST)
library(doParallel)
library(ranger)
library(sp)
library(sf)
library(terra)
library(meteo)

# prepare data
demo(meuse, echo=FALSE)
meuse <- meuse[complete.cases(meuse@data),]
data = st_as_sf(meuse, coords = c("x", "y"), crs = 28992, agr = "constant")
# data = terra::vect(meuse)
# data.frame
# data <- as.data.frame(meuse)
# data$id = 1:nrow(data)
# data.staid.x.y.z <- c(15,"x","y",NA)
fm.RFSI <- as.formula("zinc ~ dist + soil + ffreq")

# making tgrid
n.obs <- 1:6
min.node.size <- 2:10
sample.fraction <- seq(1, 0.632, -0.05) # 0.632 without / 1 with replacement
splitrule <- "variance"
ntree <- 250 # 500
mtry <- 3:(2+2*max(n.obs))
tgrid = expand.grid(min.node.size=min.node.size, num.trees=ntree,
                    mtry=mtry, n.obs=n.obs, sample.fraction=sample.fraction)

\dontrun{
# Tune RFSI model
rfsi_tuned <- tune.rfsi(formula = fm.RFSI,
                        data = data,
                        # data.staid.x.y.z = data.staid.x.y.z, # data.frame
                        # s.crs = st_crs(data),
                        # p.crs = st_crs(data),
                        tgrid = tgrid, # combinations for tuning
                        tgrid.n = 20, # number of randomly selected combinations from tgrid
                        tune.type = "LLO", # Leave-Location-Out CV
                        k = 5, # number of folds
                        seed = 42,
                        acc.metric = "RMSE", # R2, CCC, MAE
                        fit.final.model = TRUE,
                        cpus = detectCores()-1,
                        progress = TRUE,
                        importance = "impurity") # ranger parameter

rfsi_tuned$combinations
rfsi_tuned$tuned.parameters
# min.node.size num.trees mtry n.obs sample.fraction     RMSE
# 3701             3       250    6     5            0.75 222.6752
rfsi_tuned$final.model
# OOB prediction error (MSE):       46666.51 
# R squared (OOB):                  0.6517336 
}
}

