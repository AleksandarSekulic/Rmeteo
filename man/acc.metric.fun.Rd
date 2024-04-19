\name{acc.metric.fun}

\alias{acc.metric.fun}

\title{Accuracy metrics calculation}

\description{Calculates classification and regression accuracy metrics for given coresponding observation and prediction vectors.}

\usage{
acc.metric.fun(obs, pred, acc.m)
}

\arguments{
  \item{obs}{\code{numeric or factor vector;} Observations.}
  \item{pred}{\code{numeric or factor vector;} Predictions.}
  \item{acc.m}{\code{character}; Accuracy metric. Possible values for regression: "ME", "MAE", "NMAE", "RMSE", "NRMSE", "R2", "CCC". Possible values for classification: "Accuracy", "Kappa", "AccuracyLower", "AccuracyUpper", "AccuracyNull", "AccuracyPValue", "McnemarPValue".}
}

\value{
Accuracy metric value.
}

\author{Aleksandar Sekulic \email{asekulic@grf.bg.ac.rs}}

\references{
Sekulić, A., Kilibarda, M., Heuvelink, G. B., Nikolić, M. & Bajat, B. Random Forest Spatial Interpolation.Remote. Sens. 12, 1687, https://doi.org/10.3390/rs12101687 (2020).
}

\seealso{
 \code{\link{acc.metric.fun}}
 \code{\link{rfsi}}
 \code{\link{pred.rfsi}}
 \code{\link{tune.rfsi}}
 \code{\link{cv.rfsi}}
 \code{\link{pred.strk}}
 \code{\link{cv.strk}}
}

\examples{
library(sp)
library(sf)
library(CAST)
library(ranger)
library(plyr)
library(meteo)

# preparing data
demo(meuse, echo=FALSE)
meuse <- meuse[complete.cases(meuse@data),]
data = st_as_sf(meuse, coords = c("x", "y"), crs = 28992, agr = "constant")
fm.RFSI <- as.formula("zinc ~ dist + soil + ffreq")

# making tgrid
n.obs <- 1:3
min.node.size <- 2:10
sample.fraction <- seq(1, 0.632, -0.05) # 0.632 without / 1 with replacement
splitrule <- "variance"
ntree <- 250 # 500
mtry <- 3:(2+2*max(n.obs))
tgrid = expand.grid(min.node.size=min.node.size, num.trees=ntree,
                    mtry=mtry, n.obs=n.obs, sample.fraction=sample.fraction)

\donttest{
# do cross-validation
rfsi_cv <- cv.rfsi(formula=fm.RFSI, # without nearest obs
                   data = data,
                   zero.tol=0,
                   tgrid = tgrid, # combinations for tuning
                   tgrid.n = 5, # number of randomly selected combinations from tgrid for tuning
                   tune.type = "LLO", # Leave-Location-Out CV
                   k = 5, # number of folds
                   seed = 42,
                   acc.metric = "RMSE", # R2, CCC, MAE
                   output.format = "data.frame",
                   cpus=2, # detectCores()-1,
                   progress=1,
                   importance = "impurity")
summary(rfsi_cv)

# accuracy metric calculation
acc.metric.fun(rfsi_cv$obs, rfsi_cv$pred, "R2")
acc.metric.fun(rfsi_cv$obs, rfsi_cv$pred, "RMSE")
acc.metric.fun(rfsi_cv$obs, rfsi_cv$pred, "NRMSE")
acc.metric.fun(rfsi_cv$obs, rfsi_cv$pred, "MAE")
acc.metric.fun(rfsi_cv$obs, rfsi_cv$pred, "NMAE")
acc.metric.fun(rfsi_cv$obs, rfsi_cv$pred, "CCC")
}
}

