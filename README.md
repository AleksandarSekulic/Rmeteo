# `meteo`: RFSI and spacetime geostatistical interpolations for meteorological and other enviromental varialbles
[![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)
[![CRAN](http://www.r-pkg.org/badges/version/meteo)](https://cran.r-project.org/package=meteo)
[![cran checks](https://cranchecks.info/badges/worst/meteo)](https://cran.r-project.org/web/checks/check_results_meteo.html)
[![Downloads](http://cranlogs.r-pkg.org/badges/meteo?color=brightgreen)](http://www.r-pkg.org/pkg/meteo)

## Overview
`meteo` is an R package for RFSI ang spacetime geostatistical interpolations for meteorological and other enviromental varialbles.
Main functions:
* `pred.strk` - spatio-temporal regression kriging prediction (Kilibarda et al. 2014)
* `cv.strk` - k-fold cross-validation for spatio-temporal regression kriging (Kilibarda et al. 2014)
* `rfsi` - Random Forest Spatial Interpolation (RFSI) model (Sekulić et al. 2020b)
* `tune.rfsi` - tuning of RFSI model (Sekulić et al. 2020b)
* `pred.rfsi` - RFSI prediction (Sekulić et al. 2020b)
* `cv.rfsi` - nested k-fold cross-validation for RFSI (Sekulić et al. 2020b)

***Note that Out-of-bag (OOB) error statistics from RFSI model are biased and should not be considered as accuracy metrics (they do not show spatial accuracy)! The proper way to assess accuaracy of the RFSI model is by using the nested k-fold cross-validation (`cv.rfsi` function, Sekulić et al. 2020b).***

## Repositories
* [Github](https://github.com/AleksandarSekulic/Rmeteo)
* [R-forge](http://meteo.r-forge.r-project.org/)
* [CRAN](https://cran.r-project.org/package=meteo)

*Note that the latest version is in the Github repository. The R-forge and CRAN repository will updated only with a stabile version.*

## Installation
Install development versions (the most recent version) from Github with
```
library(devtools)
install_github("https://github.com/AleksandarSekulic/Rmeteo")
```
or from R-forge
```
install.packages("meteo", repos="http://R-Forge.R-project.org")
```

## Examples
### RFSI example
Complete RFSI examples (including tune.rfsi and cv.rfsi) can be found in the [demo](demo) folder.
```
library(meteo)
library(sp)
library(spacetime)
library(gstat)
library(plyr)
library(xts)
library(snowfall)
library(doParallel)
library(CAST)
library(ranger)

# preparing data
demo(meuse, echo=FALSE)
meuse <- meuse[complete.cases(meuse@data),]

#################### rfsi ####################

fm.RFSI <- as.formula("zinc ~ dist + soil + ffreq")

rfsi_model <- rfsi(formula = fm.RFSI,
                   data = meuse,
                   zero.tol = 0,
                   n.obs = 5, # number of nearest observations
                   s.crs = NA, # or meuse@proj4string # nedded only if in lon/lat (WGS84)
                   t.crs = NA, # or meuse@proj4string # nedded only if in lon/lat (WGS84)
                   cpus = detectCores()-1,
                   progress = TRUE,
                   # ranger parameters
                   importance = "impurity",
                   seed = 42,
                   num.trees = 250,
                   mtry = 5,
                   splitrule = "variance",
                   min.node.size = 5,
                   sample.fraction = 0.95,
                   quantreg = FALSE)
rfsi_model

# Note that OOB error statistics are biased and should not be considered as accuracy metrics
# (they do not show spatial accuracy)!
# The proper way to assess accuaracy of the RFSI model is by using the nested k-fold
# cross-validation (cv.rfsi function)

sort(rfsi_model$variable.importance)

#################### pred.rfsi ####################

rfsi_prediction <- pred.rfsi(model = rfsi_model,
                             data = meuse, # meuse.df (use data.staid.x.y.time)
                             zcol = "zinc",
                             newdata = meuse.grid, # meuse.grid.df (use newdata.staid.x.y.time)
                             output.format = "SpatialPixelsDataFrame",
                             zero.tol = 0,
                             s.crs = meuse@proj4string, # NA
                             newdata.s.crs = meuse@proj4string, # NA
                             t.crs = meuse@proj4string, # NA
                             cpus = detectCores()-1,
                             progress = TRUE
)

spplot(rfsi_prediction['pred'])
```
### STRK example:
Complete STRK examples (including strk.cv) can be found in the [demo](demo) folder.
```
library(meteo)
library(sp)
library(spacetime)
library(gstat)
library(plyr)
library(xts)
library(snowfall)
library(doParallel)
library(CAST)

# preparing data
data(dtempc) # temperature data
data(stations) # station locations
data(regdata) # covariates
regdata@sp@proj4string <- CRS('+proj=longlat +datum=WGS84')
data(tvgms) # ST variogram models
data(tregcoef) # MLR coefficients

serbia= point.in.polygon(stations$lon, stations$lat, c(18,22.5,22.5,18), c(40,40,46,46))
st= stations[ serbia!=0, ]
dtempc <- dtempc[dtempc$staid %in% st$staid, ]
dtempc <- dtempc[complete.cases(dtempc),]
summary(dtempc)
# create STFDF
stfdf <- meteo2STFDF(dtempc,st)
stfdf@sp@proj4string <- CRS('+proj=longlat +datum=WGS84')

#################### pred.strk ####################

# Calculate prediction of mean temperatures for "2011-07-05" and "2011-07-06" 
# global model is used for regression and variogram

results <- pred.strk(data = stfdf, # observations
                     newdata = regdata, # prediction locations with covariates
                     # newdata = regdata[,2,drop=FALSE], # for one day only
                     output.format = "stfdf",
                     reg.coef = tregcoef[[1]], # MLR coefficients
                     vgm.model = tvgms[[1]], # STRK variogram model
                     sp.nmax = 20,
                     time.nmax = 2,
                     computeVar=TRUE
)

# plot predictions
stplot(results[,,"pred", drop=F], col.regions=bpy.colors())
stplot(results[,,"var", drop=F], col.regions=bpy.colors())
```

## Citation

* Kilibarda, M., T. Hengl, G. B. M. Heuvelink, B. Graeler, E. Pebesma, M. Percec Tadic, and B. Bajat (2014), Spatio-temporal interpolation of daily temperatures for global land areas at 1 km resolution, J. Geophys. Res. Atmos., 119, 2294-2313, https://doi.org/10.1002/2013JD020803.
* Kilibarda, M., Tadić, M. P., Hengl, T., Luković, J., & Bajat, B. (2015). Global geographic and feature space coverage of temperature data in the context of spatio-temporal interpolation. Spatial Statistics, 14, 22–38. https://doi.org/10.1016/j.spasta.2015.04.005
* Sekulić, A., Kilibarda, M., Protić, D., Tadić, M. P., & Bajat, B. (2020a). Spatio-temporal regression kriging model of mean daily temperature for Croatia. Theoretical and Applied Climatology, 140(1–2), 101–114. https://doi.org/10.1007/s00704-019-03077-3
* Sekulić, A., Kilibarda, M., Heuvelink, G. B., Nikolić, M. & Bajat, B. Random Forest Spatial Interpolation.Remote. Sens. 12, 1687, https://doi.org/10.3390/rs12101687 (2020b).

