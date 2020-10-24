cv.strk <- function (data, # data.frame(id,x,y,time,obs,ec1,ec2,...) | STFDF - with covariates
                     zcol=1,
                     data.staid.x.y.time = c(1,2,3,4), # if data.frame
                     obs, # data.frame(id,time,obs,cov)
                     obs.staid.time = c(1,2),
                     stations, # data.frame(id,x,y)
                     stations.staid.x.y = c(1,2,3),
                     zero.tol=0,
                     reg.coef, # check coef names
                     vgm.model,
                     sp.nmax=20, # use all if not specified
                     time.nmax=2, # use all if not specified
                     type = "LLO", # type of cv - LLO for now, after LTO, LLTO - CAST
                     k = 5, # number of folds
                     seed = 42,
                     folds, # if user wants to create folds
                     fold.column, # by which column
                     refit = TRUE,
                     output.format = "STFDF",
                     parallel.processing = FALSE, # doParallel
                     pp.type = "snowfall", # "doParallel"
                     cpus=detectCores()-1,
                     progress=TRUE,
                     ...){
  
  # check the input
  print('Preparing data ...')
  if ((missing(data) & missing(obs) & missing(stations)) | missing(reg.coef) | missing(vgm.model)) {
    stop('The arguments data (or obs and stations), reg.coef and vgm.model must not be empty!')
  }
  if (!missing(data)){
    if (class(data) == "data.frame") {
      # if zcol is character
      if (is.numeric(zcol)) {
        zcol.name <- names(data)[zcol]
      } else {
        zcol.name <- zcol
        zcol <- index(names(data))[names(data) == zcol.name]
      }
      # if data.staid.x.y.time is character
      if (!is.numeric(data.staid.x.y.time)) {
        data.staid.x.y.time <- sapply(data.staid.x.y.time, function(i) index(names(data))[names(data) == i])
      }
      # to stfdf
      obs <- cbind(data[, c(data.staid.x.y.time[c(1,4)], zcol)], data[, -c(data.staid.x.y.time[c(1,4)], zcol)])
      if ('endTime' %in% names(obs)) { obs <- obs[, -which('endTime' == names(obs))] }
      if ('timeIndex' %in% names(obs)) { obs <- obs[, -which('timeIndex' == names(obs))] }
      stations <- data[, c(data.staid.x.y.time[1:3])]
      stations <- unique(stations[complete.cases(stations), ])
      data <- meteo2STFDF(obs      = obs,
                          stations = stations,
                          crs      = CRS("+proj=longlat +datum=WGS84"),
                          obs.staid.time = c(1,2),
                          stations.staid.lon.lat = c(1,2,3)
      )
      zcol=1
    } else if (class(data) == "STFDF" | class(data) == "STSDF") { # | class(data) != "STSDF
      if (is.numeric(zcol)) {
        zcol.name <- names(data@data)[zcol]
      } else {
        zcol.name <- zcol
        zcol <- index(names(data@data))[names(data@data) == zcol.name]
      }
    } else {
      stop('The argument data must be of STFDF, STSDF or data.frame class!') # "STSDF"
    }
  } else { # obs, stations
    # if zcol is character
    if (is.numeric(zcol)) {
      zcol.name <- names(obs)[zcol]
    } else {
      zcol.name <- zcol
      zcol <- index(names(obs))[names(obs) == zcol.name]
    }
    # if data.staid.x.y.time is character
    if (!is.numeric(obs.staid.time)) {
      obs.staid.time <- sapply(obs.staid.time, function(i) index(names(obs))[names(obs) == i])
    }
    # if stations.staid.x.y is character
    if (!is.numeric(stations.staid.x.y)) {
      stations.staid.x.y <- sapply(stations.staid.x.y, function(i) index(names(stations))[names(stations) == i])
    }
    # to stfdf
    obs <- cbind(obs[, c(obs.staid.time, zcol)], obs[, -c(obs.staid.time, zcol)])
    if ('endTime' %in% names(obs)) { obs <- obs[, -which('endTime' == names(obs))] }
    if ('timeIndex' %in% names(obs)) { obs <- obs[, -which('timeIndex' == names(obs))] }
    stations <- stations[, c(stations.staid.x.y)]
    stations <- unique(stations[complete.cases(stations), ])
    data <- meteo2STFDF(obs      = obs,
                        stations = stations,
                        crs      = CRS("+proj=longlat +datum=WGS84"),
                        obs.staid.time = c(1,2),
                        stations.staid.lon.lat = c(1,2,3)
    )
    zcol=1
  }

  time <- index(data@time)
  data.df <- as.data.frame(data)
  names_covar <- names(reg.coef)[-1]
  
  c.dif <- setdiff(names_covar, names(data.df))
  
  # if covariate names doesn't exist in data - do ordinary kriging
  if (!identical(c.dif, character(0))) {
    # if covariates - do overlay
    warning(paste('The covariate(s) ', paste(c.dif, collapse = ", "), ' - missing from data!', sep = ""))
    warning('Trend is set to 0, performing space-time ordinary kriging cross-validation.')
    # data$tlm<-0 
    # data$tres <- data@data[,zcol]
    names_covar <- names(data@data[zcol])[1]
  }
  
  # set folds
  if (missing(fold.column)){
    if (missing(folds)) {
      # create folds
      space_id <- rep(1:length(data@sp), length(time))
      time_id <- rep(1:length(time), each = length(data@sp))
      st_df <- cbind(space_id, time_id)
      # if (type == "LLO") {
        spacevar <- "space_id"
        timevar <- NA
      # TO DO LTO and LLTO
      # } else if (type == "LTO") {
      #   spacevar <- NA
      #   timevar <- "time_id"
      # } else if (type == "LLTO") {
      #   spacevar <- "space_id"
      #   timevar <- "time_id"
      # }
      indices <- CreateSpacetimeFolds(st_df, spacevar = spacevar, timevar = timevar,
                                      k=k, seed = seed)
      folds <- c()
      for (f in 1:length(indices$indexOut)) {
        folds[indices$indexOut[[f]]] <- f
      }
    }
    data$folds <- folds
    fold.column <- "folds"
  }
  
  # remove duplicates
  data <- rm.dupl(data, zcol, zero.tol)
  
  # remove the stations where covariate is missing
  nrowsp <- length(data@sp)
  for (covar in names_covar){
    # count NAs per stations
    if (covar %in% names(data@data)) {
      numNA <- apply(matrix(data@data[,covar],
                            nrow=nrowsp,byrow=F), MARGIN=1,
                     FUN=function(x) sum(is.na(x)))
      # Remove stations out of covariates
      rem <- numNA != length(time)
      data <-  data[rem,drop=F]
    } else {
      data@sp <- data@sp[!is.na(data@sp@data[, covar]), ]
    }
  }
  
  # Remove dates out of covariates
  rm.days <- c()
  for (t in 1:length(time)) {
    if(sum(complete.cases(data[, t]@data)) == 0) {
      rm.days <- c(rm.days, t)
    }
  }
  if(!is.null(rm.days)){
    data <- data[,-rm.days]
  }
  
  data.df <- as.data.frame(data)
  
  if(nrow(data.df[complete.cases(data.df), ]) == 0){
    stop('The data does not have complete cases!')
  } else {
    print(paste(nrow(data.df[complete.cases(data.df), ]), " observations complete cases (", length(data@sp), " stations x ", length(data@time), " days) used for cross-validation.", sep=""))
  }
  
  # do CV
  print('Doing CV ...')
  pred <- c()
  vgm.model.init <- vgm.model
  for (val_fold in sort(unique(data.df[, fold.column]))) {
    print(paste('Fold ', val_fold, sep=""))
    
    dev <- data
    dev@data[dev$folds==val_fold & !is.na(dev$folds), ] <- NA
    dev <- as(dev, "STSDF")
    val <- data
    val@data[val$folds!=val_fold & !is.na(val$folds), ] <- NA
    val <- as(val, "STSDF")
    dev_df <- as.data.frame(dev)
    val_df <- as.data.frame(val)
    # dev_df = data.df[data.df$folds!=val_fold, ]
    # val_df = data.df[data.df$folds==val_fold, ]
    
    fcv <- paste(zcol.name, " ~ ", paste(names_covar, collapse = " + "), sep="")
    dev_df <- dev_df[, c(zcol.name, names_covar)]
    dev_df$completed = complete.cases(dev_df)
    
    if(refit) {
      # fit lm - use variables in reg.coef + update reg.coef
      ### fiting STRK model on dev data ###
      set.seed(seed)
      fold_lm = lm(fcv, dev_df[dev_df$completed, ])
      reg.coef = coefficients(fold_lm) # model coefficients, added to tregcoef$tmeanHR
      # summary(fold_lm)
      lm_trend = c()
      br = 1
      for (i in 1:length(dev_df$completed)) {
        if (dev_df$completed[i]){
          lm_trend[i] = fold_lm$fitted.values[br]
          br = br + 1
        } else {
          lm_trend[i] = NA
        }
      }
      dev$lm_trend = lm_trend
      dev$lm_res = dev@data[, zcol] - dev$lm_trend
      
      # rescale varigram distance
      var = variogramST(lm_res ~ 1, dev) #, ...) # , tlags = 0:5, cutoff = 300, width = 10, na.omit=T) # tunit="days"
      
      scale.num <- round(max(var$dist/10, na.rm=T) / 10) * 10
      if (scale.num == 0) {
        scale.num = 1
      }
      var$dist = var$dist/scale.num
      var$spacelag = var$spacelag/scale.num
      var$avgDist = var$avgDist/scale.num
      attr(var, "boundaries") = attr(var, "boundaries") / scale.num
      
      vgm.model <- vgm.model.init
      vgm.model$space$range <- vgm.model$space$range / scale.num
      vgm.model$joint$range <- vgm.model$joint$range / scale.num
      vgm.model$stAni <- vgm.model$stAni / scale.num
      
      pars.l <- c(sill.s = 0, range.s = 1, nugget.s = 0,
                  sill.t = 0, range.t = 1, nugget.t = 0,
                  sill.st = 0, range.st = 1, nugget.st = 0,
                  anis = 0)
      
      # fit variogram - use vgm.model as initial + update vgm.model
      sumMetric_Vgm <- fit.StVariogram(var, vgm.model, method="L-BFGS-B",lower=pars.l)
      # attr(sumMetric_Vgm, "MSE")
      # sumMetric_Vgm
      
      sumMetric_Vgm$space$range <- sumMetric_Vgm$space$range * scale.num
      sumMetric_Vgm$joint$range <- sumMetric_Vgm$joint$range * scale.num
      sumMetric_Vgm$stAni <- sumMetric_Vgm$stAni * scale.num
      vgm.model <- sumMetric_Vgm
      
    }
    
    # validation #
    val_stfdf <- pred.strk(data = dev, # stfdf.df
                           zcol = zcol.name,
                           newdata = val, # regdata.df
                           output.format = "STSDF", # df
                           zero.tol = zero.tol,
                           reg.coef = reg.coef, # check coef names
                           vgm.model = vgm.model,
                           sp.nmax = sp.nmax, # use all if not specified
                           time.nmax = time.nmax, # use all if not specified
                           by = 'time',
                           parallel.processing = parallel.processing, # doParallel
                           pp.type = pp.type, # "doParallel"
                           cpus = cpus,
                           computeVar = F,
                           ...
    )
    
    val_df <- as.data.frame(val)[, c(names(val_df)[1:4], fold.column, zcol.name)]
    val_df$pred <- val_stfdf$pred
    pred <- rbind(pred, val_df)
  }
  names(pred)[names(pred) == zcol.name] <- "obs"
  
  # return
  if (output.format == "data.frame") {
    print("Done!")
    return(pred)
  } else {
    sta <- pred[, c(3,1,2)]
    sta <- sta[!duplicated(sta), ]
    obs <- pred[, 3:length(pred)]
    stfdf <- meteo2STFDF(obs      = obs,
                         stations = sta,
                         crs      = CRS("+proj=longlat +datum=WGS84"),
                         obs.staid.time = c(1,2),
                         stations.staid.lon.lat = c(1,2,3)
    )
    if (output.format == "STSDF") {
    print("Done!")
    return(as(stfdf, "STSDF"))
    } else { #  (output.format == "STFDF")
      print("Done!")
      return(stfdf)
    }
  }
  
  ##########
  
}
