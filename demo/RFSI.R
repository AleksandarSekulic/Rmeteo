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

#################### near.obs ####################
nearest_obs <- near.obs(locations = meuse, # from
                        locations.x.y = c("x","y"),
                        observations = meuse, # to
                        observations.x.y=c("x","y"),
                        zcol = "zinc",
                        n.obs = 10, # number of nearest observations
                        rm.dupl = TRUE) # uif TRUE, the observations themself will be excluded from spatial covariates 
str(nearest_obs)
summary(nearest_obs)

#################### rfsi ####################

# fm.RFSI <- as.formula("om ~ dist + soil") # covariates without spatial covariates (obs1, dist1, ...)
fm.RFSI <- as.formula("zinc ~ dist + soil + ffreq")
meuse.df <- as.data.frame(meuse)
meuse.df$id = 1:nrow(meuse.df)

rfsi_model <- rfsi(formula = fm.RFSI,
                   data = meuse, # meuse.df (use data.staid.x.y.time)
                   # data.staid.x.y.time = c("id", "x", "y", NA), # only if class(data) == data.frame
                   # obs = meuse.df, # uses obs and stations in combination
                   # obs.staid.time = c("id", NA),
                   # stations = meuse.df,
                   # stations.staid.x.y = c("id", "x", "y"),
                   zero.tol = 0,
                   n.obs = 5, # number of nearest observations
                   s.crs = NA, # or meuse@proj4string # nedded only if the coordinates are lon/lat (WGS84)
                   t.crs = NA, # or meuse@proj4string # nedded only if the coordinates are lon/lat (WGS84)
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
# Note that OOB error statistics are biased and should not be considered as accuracy metrics (they do not show spatial accuracy)!
# The proper way to assess accuaracy of the RFSI model is by using the nested k-fold cross-validation (cv.rfsi function)
sort(rfsi_model$variable.importance)

#################### pred.rfsi ####################

meuse.grid.df <- as.data.frame(meuse.grid)
meuse.grid.df$id <- 1:nrow(meuse.grid.df)

rfsi_prediction <- pred.rfsi(model = rfsi_model,
                             data = meuse, # meuse.df (use data.staid.x.y.time)
                             zcol = "zinc",
                             # data.staid.x.y.time = c("id", "x", "y", NA), # only if class(data) == data.frame
                             # obs = meuse.grid,
                             # obs.staid.time = c("id", NA),
                             # stations = meuse.grid,
                             # stations.staid.x.y = c("id", "x", "y"),
                             newdata = meuse.grid, # meuse.grid.df (use newdata.staid.x.y.time)
                             # newdata.staid.x.y.time = c("id", "x", "y", NA), # only if class(newdata) == data.frame
                             output.format = "SpatialPixelsDataFrame",
                             zero.tol = 0,
                             s.crs = meuse@proj4string, # NA
                             newdata.s.crs = meuse@proj4string, # NA
                             t.crs = meuse@proj4string, # NA
                             cpus = detectCores()-1,
                             progress = TRUE
)
spplot(rfsi_prediction['pred'])

#################### tune.rfsi ####################

# making tgrid
n.obs <- 2:6
min.node.size <- 2:10
sample.fraction <- seq(1, 0.632, -0.05) # 0.632 without / 1 with replacement
splitrule <- "variance"
ntree <- 250 # 500
mtry <- 3:(2+2*max(n.obs))
tgrid = expand.grid(min.node.size=min.node.size, num.trees=ntree,
                    mtry=mtry, n.obs=n.obs, sample.fraction=sample.fraction)

# Tune an RFSI model
rfsi_tuned <- tune.rfsi(formula = fm.RFSI,
                        data = meuse, # meuse.df (use data.staid.x.y.time)
                        # data.staid.x.y.time = c("id", "x", "y", NA), # only if class(data) == data.frame
                        zero.tol = 0,
                        s.crs = NA,
                        t.crs = NA,
                        tgrid = tgrid, # combinations for tuning
                        tgrid.n = 5, # number of randomly selected combinations from tgrid for tuning
                        tune.type = "LLO", # Leave-Location-Out CV
                        k = 5, # number of folds
                        acc.metric = "RMSE", # R2, CCC, MAE
                        fit.final.model = TRUE,
                        cpus = detectCores()-1,
                        progress = TRUE,
                        # ranger parameters
                        importance = "impurity",
                        seed = 42)

rfsi_tuned$combinations
rfsi_tuned$tuned.parameters
rfsi_tuned$final.model

#################### cv.rfsi ####################

# Cross-validation of RFSI
rfsi_cv <- cv.rfsi(formula=fm.RFSI, # without nearest obs
                   data = meuse, # meuse.df (use data.staid.x.y.time)
                   # data.staid.x.y.time = c("id", "x", "y", NA), # only if class(data) == data.frame
                   zero.tol=0,
                   s.crs=NA,
                   t.crs=NA,
                   tgrid = tgrid, # combinations for tuning
                   tgrid.n = 5, # number of randomly selected combinations from tgrid for tuning
                   tune.type = "LLO", # Leave-Location-Out CV
                   k = 5, # number of folds
                   seed = 42,
                   acc.metric = "RMSE", # R2, CCC, MAE
                   output.format = "data.frame", # "SpatialPointsDataFrame",
                   cpus=detectCores()-1,
                   progress=FALSE,
                   # ranger parameters
                   importance = "impurity")

summary(rfsi_cv)
# spplot(rfsi_cv[, , "pred"])
# spplot(rfsi_cv[, , "obs"])
acc.metric.fun(rfsi_cv$obs, rfsi_cv$pred, "R2")
acc.metric.fun(rfsi_cv$obs, rfsi_cv$pred, "RMSE")
acc.metric.fun(rfsi_cv$obs, rfsi_cv$pred, "MAE")
acc.metric.fun(rfsi_cv$obs, rfsi_cv$pred, "CCC")
