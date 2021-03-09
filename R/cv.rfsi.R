cv.rfsi <- function (formula, # without nearest obs
                     data, # data.frame(x,y,obs,time,ec1,ec2,...) | STFDF - with covariates | SpatialPointsDataFrame | SpatialPixelsDataFrame
                     data.staid.x.y.time = c(1,2,3,4), # if data.frame
                     obs, # data.frame(id,time,obs,cov)
                     obs.staid.time = c(1,2),
                     stations, # data.frame(id,x,y)
                     stations.staid.x.y = c(1,2,3),
                     zero.tol=0,
                     use.idw=FALSE,
                     # avg = FALSE,
                     # increment, # avg(nearest point dist)
                     # range, # bbox smaller(a, b) / 2
                     # direct = FALSE,
                     s.crs=NA,
                     t.crs=NA,
                     tgrid, # caret tune grid - by default random (min.node.size, mtry, no, sample.fraction, ntree, splitrule)
                     tgrid.n=10,
                     tune.type = "LLO", # type of cv - LLO for now, after LTO, LLTO - CAST
                     k = 5, # number of folds
                     seed=42,
                     folds, # if user want to create folds
                     fold.column, # by which column
                     acc.metric, # for tuning on subfolds
                     output.format = "data.frame", #"STFDF",
                     cpus=detectCores()-1,
                     progress=TRUE,
                     soil3d = FALSE, # soil RFSI
                     no.obs = 'increase', # exactly
                     ...){ # fixed ranger parameters
  
  # check the input
  if (progress) print('Preparing data ...')
  if ((missing(data) & missing(obs) & missing(stations)) | missing(formula) | missing(tgrid)) {
    stop('The arguments data (or obs and stations), formula and tgrid must not be empty!')
  }
  formula <- as.formula(formula)
  all_vars <- all.vars(formula)
  zcol.name <- all_vars[1]
  names_covar <- all_vars[-1]
  
  if (!missing(data)){
    if (class(data) == "data.frame") {
      # if data.staid.x.y.time is character
      if (!is.numeric(data.staid.x.y.time)) {
        data.staid.x.y.time <- match(data.staid.x.y.time, names(data)) #sapply(data.staid.x.y.time, function(i) index(names(data))[names(data) == i])
      }
      data.df = data
    } else if (class(data) == "STFDF" | class(data) == "STSDF") {
      if (class(data) == "STSDF") {data <- as(data, "STFDF")}
      data <- rm.dupl(data, zcol, zero.tol)
      data.df <- as.data.frame(data)
      data.staid.x.y.time <- c(3,1,2,4)
      if (!is.na(data@sp@proj4string)) {
        s.crs <- data@sp@proj4string
      }
    } else if (class(data) == "SpatialPointsDataFrame" | class(data) == "SpatialPixelsDataFrame") {
      data.df <- as.data.frame(data)
      data.df$staid <- 1:nrow(data.df)
      data.staid.x.y.time <- c(length(data.df),length(data.df)-2,length(data.df)-1,NA)
      if (!is.na(data@sp@proj4string)) {
        s.crs <- data@proj4string
      }
    } else {
      stop('The argument data must be of STFDF, STSDF, data.frame,SpatialPointsDataFrame or SpatialPixelsDataFrame class!') # "STSDF"
    }
  } else { # obs, stations
    # if obs.staid.time is character
    if (!is.numeric(obs.staid.time)) {
      obs.staid.time <- match(obs.staid.time, names(obs)) # sapply(obs.staid.time, function(i) index(names(obs))[names(obs) == i])
    }
    # if stations.staid.x.y is character
    if (!is.numeric(stations.staid.x.y)) {
      stations.staid.x.y <- match(stations.staid.x.y, names(stations)) # sapply(stations.staid.x.y, function(i) index(names(stations))[names(stations) == i])
    }
    # to stfdf
    data.df <- join(obs, stations, by=names(obs)[obs.staid.time[1]], match="first")
    data.staid.x.y.time <- c(obs.staid.time[1],
                             stations.staid.x.y[2] + length(obs),
                             stations.staid.x.y[3] + length(obs),
                             obs.staid.time[2])
  }
  
  # Criteria accuracy parameter
  if (is.factor(data.df[, zcol.name])) {
    # classification
    if (missing(acc.metric) || !(acc.metric %in% c("Accuracy","Kappa","AccuracyLower","AccuracyUpper","AccuracyNull","AccuracyPValue","McnemarPValue"))) {
      acc.metric <- "Kappa"
      warning("Criteria accuracy parameter is missing or is not valid for classification task. Using Kappa acccuracy parameter!")
    }
  } else {
    # regression
    if (missing(acc.metric) || !(acc.metric %in% c("RMSE","MAE","ME","R2","CCC"))) {
      acc.metric <- "RMSE"
      warning("Criteria accuracy parameter is missing or is not valid for regression task. Using RMSE acccuracy parameter!")
    }
  }
  
  # set folds
  if (missing(fold.column)){
    if (missing(folds)) {
      # create folds
      
      # check if time?
      
      # space_id <- rep(1:length(data@sp), length(time))
      # time_id <- rep(1:length(time), each = length(data@sp))
      # st_df <- cbind(space_id, time_id)
      # if (type == "LLO") {
      spacevar <- names(data.df)[data.staid.x.y.time[1]]# "space_id"
      timevar <- NA
      # TO DO LTO and LLTO
      # } else if (type == "LTO") {
      #   spacevar <- NA
      #   timevar <- "time_id"
      # } else if (type == "LLTO") {
      #   spacevar <- "space_id"
      #   timevar <- "time_id"
      # }
      indices <- CreateSpacetimeFolds(data.df, spacevar = spacevar, timevar = timevar,
                                      k=k, seed = seed)
      folds <- c()
      for (f in 1:length(indices$indexOut)) {
        folds[indices$indexOut[[f]]] <- f
      }
    }
    data.df$folds <- folds
    fold.column <- "folds"
  }
  
  # do CV
  if (progress) print('Doing CV ...')
  pred <- c()
  # for val_fold in main folds
  for (val_fold in sort(unique(data.df[, fold.column]))) {
    print(paste('### Main fold ', val_fold, " ###", sep=""))
    # tune RFSI model
    dev.df <- data.df[data.df[, fold.column] != val_fold, ]
    val.df <- data.df[data.df[, fold.column] == val_fold, ]
    if (progress) print('Tuning RFSI model ...')
    # tune.rfsi
    tuned_model <- tune.rfsi(formula, # without nearest obs
                             data=dev.df, # data.frame(x,y,obs,time,ec1,ec2,...) | STFDF - with covariates | SpatialPointsDataFrame | SpatialPixelsDataFrame
                             data.staid.x.y.time = data.staid.x.y.time, # if data.frame
                             zero.tol=zero.tol,
                             # n.obs=n.obs, # nearest obs
                             # time.nmax, # use all if not specified
                             s.crs=s.crs,
                             t.crs=t.crs,
                             # use.tps=use.tps,
                             # tps.df=tps.df,
                             use.idw=use.idw,
                             # parallel.processing = FALSE, # doParallel - videti zbog ranger-a
                             tgrid=tgrid, # caret tune grid - by default random (min.node.size, mtry, no, sample.fraction, ntree, splitrule)
                             tgrid.n=tgrid.n,
                             tune.type = tune.type, # type of cv - LLO for now, after LTO, LLTO - CAST
                             k = k, # number of folds
                             # folds=folds, # if user wants to create folds
                             # fold.column=fold.column, # by which column
                             acc.metric = acc.metric,
                             cpus=cpus, # for near.obs
                             progress=progress,
                             fit.final.model=TRUE,
                             soil3d = soil3d, # soil RFSI
                             no.obs = no.obs, # exactly
                             # quantreg = quantreg,
                             # seed = seed,
                             ...)
    # validate
    if (progress) print('Validation ...')
    fold_prediction <- pred.rfsi(tuned_model$final.model, # RFSI model iz rfsi ili tune rfsi funkcije
                                 data=dev.df, # data.frame(x,y,obs,time) | STFDF - with covariates | SpatialPointsDataFrame | SpatialPixelsDataFrame
                                 zcol=zcol.name,
                                 data.staid.x.y.time = data.staid.x.y.time, # if data.frame
                                 newdata=val.df, # data.frame(x,y,time,ec1,ec2,...) | STFDF - with covariates | SpatialPointsDataFrame | SpatialPixelsDataFrame
                                 newdata.staid.x.y.time = data.staid.x.y.time, # if data.frame
                                 output.format = "data.frame",
                                 zero.tol=zero.tol,
                                 # n.obs=10, # nearest obs 3 vidi iz modela
                                 # time.nmax, # use all if not specified
                                 s.crs=s.crs,
                                 newdata.s.crs = s.crs,
                                 t.crs=t.crs,
                                 # parallel.processing = FALSE, # doParallel - videti zbog ranger-a
                                 cpus=cpus,
                                 progress=progress,
                                 soil3d = soil3d, # soil RFSI
                                 depth.range = tuned_model$tuned.parameters$depth.range, # in units of depth
                                 no.obs = no.obs, # exactly
                                 ...)
    
    if (is.na(data.staid.x.y.time[4])) {
      fold_prediction <- fold_prediction[, c(names(val.df)[data.staid.x.y.time[1:3]], "pred")]
      val.df <- val.df[, c(names(val.df)[data.staid.x.y.time[1:3]], fold.column, zcol.name)]
    } else {
      fold_prediction <- fold_prediction[, c(names(val.df)[data.staid.x.y.time], "pred")]
      val.df <- val.df[, c(names(val.df)[data.staid.x.y.time], fold.column, zcol.name)]
    }
    ### ovde radi join kao kod tune!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    val.df <- join(val.df, fold_prediction)
    # val.df$pred <- fold_prediction$pred # ovo proveri za quantreg???
    pred <- rbind(pred, val.df)
  
    # fold_obs <- c(fold_obs, val.df[, zcol.name])
    # fold_pred <- c(fold_pred, fold_prediction$pred)
  }
  names(pred)[names(pred) == zcol.name] <- "obs"
  pred <- pred[complete.cases(pred), ]
  
  # return
  if (output.format == "STSDF" | output.format == "STFDF" ) {
    sta <- pred[, 1:3]
    sta <- sta[!duplicated(sta), ]
    obs <- pred[, 4:length(pred)]
    stfdf <- meteo2STFDF(obs      = obs,
                         stations = sta,
                         crs      = s.crs,
                         obs.staid.time = c(1,2),
                         stations.staid.lon.lat = c(1,2,3)
    )
    if (output.format == "STSDF") {
      if (progress) print("Done!")
      return(as(stfdf, "STSDF"))
    } else { #  (output.format == "STFDF")
      if (progress) print("Done!")
      return(stfdf)
    }
  } else if (output.format == "SpatialPointsDataFrame") {
    spdf <- pred
    coordinates(spdf) <- c(names(spdf)[2:3])
    spdf@proj4string <- s.crs
    if (progress) print("Done!")
    return(spdf)
  } else if (output.format == "SpatialPixelsDataFrame") {
    spdf <- pred
    coordinates(spdf) <- c(names(spdf)[2:3])
    spdf@proj4string <- s.crs
    spdf <- as(spdf, "SpatialPixelsDataFrame")
    if (progress) print("Done!")
    return(spdf)
  } else { #  (output.format == "df")
    if (progress) print("Done!")
    return(pred)
  }
  
}

# list - final RF model, parameters - n.obs, mtry, min.node.size, splitrule