\name{cv.rfsi}

\alias{cv.rfsi}

\title{Nested k-fold cross-validation for Random Forest Spatial Interpolation (RFSI)}

\description{Function for nested k-fold cross-validation function for Random Forest Spatial Interpolation (RFSI) (Sekulić et al. 2020). It is based on \link{rfsi}, \link{pred.rfsi}, and \link{tune.rfsi} functions. Currently, only spatial (leave-location-out) cross-validation is implemented. Temporal and spatio-temporal cross-validation will be implemented in the future.}

\usage{
cv.rfsi(formula,
        data,
        data.staid.x.y.z = NULL,
        zero.tol = 0,
        use.idw = FALSE,
        s.crs = NA,
        p.crs = NA,
        tgrid,
        tgrid.n=10,
        tune.type = "LLO",
        k = 5,
        seed=42,
        folds,
        fold.column,
        acc.metric,
        output.format = "data.frame",
        cpus = detectCores()-1,
        progress = 1,
        soil3d = FALSE,
        no.obs = 'increase',
        ...)
}

\arguments{
  \item{formula}{formula; Formula for specifying target variable and covariates (without nearest observations and distances to them). If \code{z~1}, an RFSI model using only nearest obsevrations and distances to them as covariates will be cross-validated.}
  \item{data}{\link[sf]{sf-class}, \link[sftime]{sftime-class}, \link[terra]{SpatVector-class} or \code{data.frame}; Contains target variable (observations) and covariates used for making an RFSI model. If \code{data.frame} object, it should have next columns: station ID (staid), longitude (x), latitude (y), 3rd component - time, depth, ... (z) of the observation, observation value (obs) and covariates (cov1, cov2, ...). If covariates are missing, the RFSI model using only nearest obsevrations and distances to them as covariates (\code{formula=z~1}) will be cross-validated.}
  \item{data.staid.x.y.z}{numeric or character vector; Positions or names of the station ID (staid), longitude (x), latitude (y) and 3rd component (z) columns in data.frame \code{data} object (e.g. c(1,2,3,4)). If \code{data} is \link[sf]{sf-class}, \link[sftime]{sftime-class}, or \link[terra]{SpatVector-class} object, \code{data.staid.x.y.z} is used to point staid and z position. Set z position to NA (e.g. c(1,2,3,NA)) or ommit it (e.g. c(1,2,3)) for spatial interpolation. Default is NULL.}
  \item{zero.tol}{numeric; A distance value below (or equal to) which locations are considered as duplicates. Default is 0. See \link{rm.dupl}.}
  \item{use.idw}{boolean; IDW prediction as covariate - will IDW predictions from \code{n.obs} nearest observations be calculated and tuned (see function \link{near.obs}). Default is FALSE.}
  \item{s.crs}{\link[sf]{crs-class}; Source CRS of \code{data}. If \code{data} contains crs, \code{s.crs} will be overwritten. Default is NA.}
  \item{p.crs}{\link[sf]{crs-class}; Projection CRS for \code{data} reprojection. If NA, \code{s.crs} will be used for distance calculation. Note that observations should be in projection for finding nearest observations based on Eucleadean distances (see function \link{near.obs}). Default is NA.}
  \item{tgrid}{data.frame; Possible tuning parameters for nested folds. The column names are same as the tuning parameters. Possible tuning parameters are: \code{n.obs}, \code{num.trees}, \code{mtry}, \code{min.node.size}, \code{sample.fraction}, \code{splirule}, \code{idw.p}, and \code{depth.range}.}
  \item{tgrid.n}{numeric; Number of randomly chosen \code{tgrid} combinations for nested tuning of RFSI. If larger than \code{tgrid}, will be set to \code{length(tgrid)}}
  \item{tune.type}{character; Type of nested cross-validation: leave-location-out ("LLO"), leave-time-out ("LTO") - TO DO, and leave-location-time-out ("LLTO") - TO DO. Default is "LLO".}
  \item{k}{numeric; Number of random outer and inner folds (i.e. for cross-validation and nested tuning) that will be created with \link[CAST]{CreateSpacetimeFolds} function. If \code{folds} or \code{fold.column} parameters are specified, then used only for inner folds. Default is 5.}
  \item{seed}{numeric; Random seed that will be used to generate outer  and inner folds with \link[CAST]{CreateSpacetimeFolds} function.}
  \item{folds}{numeric or character vector; Showing outer folds of \code{data} observations (rows) used for cross-validation. Used if \code{fold.column} parameter is not specified.  Note that inner folds for nested tuning will be created with \link[CAST]{CreateSpacetimeFolds} function.}
  \item{fold.column}{numeric or character; Column name or number showing the position of column in \code{data} that represents outer folds used for cross-validation. Note that inner folds for nested tuning will be created with \link[CAST]{CreateSpacetimeFolds} function.}
  \item{acc.metric}{character; Accuracy metric that will be used as a criteria for choosing an optimal RFSI model in nested tuning. Possible values for regression: "ME", "MAE", "NMAE", "RMSE" (default), "NRMSE", "R2", "CCC". Possible values for classification: "Accuracy","Kappa" (default), "AccuracyLower", "AccuracyUpper", "AccuracyNull", "AccuracyPValue", "McnemarPValue".}
  \item{output.format}{character; Format of the output, \code{data.frame} (default), \link[sf]{sf-class}, \link[sftime]{sftime-class}, or \link[terra]{SpatVector-class}.}
  \item{cpus}{numeric; Number of processing units. Default is detectCores()-1.}
  \item{progress}{numeric; If progress bar is shown. 0 is no progress bar, 1 is outer folds results, 2 is + innner folds results, 3 is + prediction progress bar. Default is 1.}
  \item{soil3d}{logical; If 3D soil modellig is performed and \link{near.obs.soil} function is used for finding n nearest observations and distances to them. In this case, z position of the \code{data.staid.x.y.z} points to the depth column.}
  \item{no.obs}{character; Possible values are \code{increase} (default) and \code{exactly}. If set to \code{increase}, in case if there is no \code{n.obs} observations in \code{depth.range} for a specific location, the \code{depth.range} is increased (multiplied by 2, 3, ...) until the number of observations are larger or equal to \code{n.obs}. If set to \code{exactly}, the function will raise an error when it come to the first location with no \code{n.obs} observations in specified \code{depth.range} (see function \link{near.obs.soil}).}
  \item{...}{Further arguments passed to \link[ranger]{ranger}.}
}

\value{
A \code{data.frame}, \link[sf]{sf-class}, \link[sftime]{sftime-class}, or \link[terra]{SpatVector-class} obejct (depends on \code{output.format} argument), with columns:
  \item{obs}{Observations.}
  \item{pred}{Predictions from cross-validation.}
  \item{folds}{Folds used for cross-validation.}
}

\author{Aleksandar Sekulic \email{asekulic@grf.bg.ac.rs}}

\references{
Sekulić, A., Kilibarda, M., Heuvelink, G. B., Nikolić, M. & Bajat, B. Random Forest Spatial Interpolation. Remote. Sens. 12, 1687, https://doi.org/10.3390/rs12101687 (2020).
}

\seealso{
  \code{\link{near.obs}}
  \code{\link{rfsi}}
  \code{\link{pred.rfsi}}
  \code{\link{tune.rfsi}}
}

\examples{
library(CAST)
library(doParallel)
library(dplyr)
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

# Cross-validation of RFSI
rfsi_cv <- cv.rfsi(formula=fm.RFSI, # without nearest obs
                   data = data,
                   # data.staid.x.y.z = c("id", "x", "y", NA), # only if class(data) == data.frame
                   zero.tol=0,
                   # s.crs=NA,
                   # t.crs=NA,
                   tgrid = tgrid, # combinations for tuning
                   tgrid.n = 5, # number of randomly selected combinations from tgrid for tuning
                   tune.type = "LLO", # Leave-Location-Out CV
                   k = 5, # number of folds
                   seed = 42,
                   acc.metric = "RMSE", # R2, CCC, MAE
                   output.format = "sf", # "data.frame", # "SpatVector",
                   cpus=detectCores()-1,
                   progress=1,
                   importance = "impurity") # ranger parameter

summary(rfsi_cv)
rfsi_cv$dif <- rfsi_cv$obs - rfsi_cv$pred
plot(rfsi_cv["dif"])
# spplot(rfsi_cv[, , "obs"])
acc.metric.fun(rfsi_cv$obs, rfsi_cv$pred, "R2")
acc.metric.fun(rfsi_cv$obs, rfsi_cv$pred, "RMSE")
acc.metric.fun(rfsi_cv$obs, rfsi_cv$pred, "MAE")
acc.metric.fun(rfsi_cv$obs, rfsi_cv$pred, "CCC")
}

