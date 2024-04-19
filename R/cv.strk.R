cv.strk <- function (data, # data.frame(id,x,y,time,obs,ec1,ec2,...) | STFDF - with covariates
                     obs.col=1,
                     data.staid.x.y.z = NULL, # if data.frame
                     crs = NA, # brisi
                     zero.tol=0,
                     reg.coef, # check coef names
                     vgm.model,
                     sp.nmax=20, # use all if not specified
                     time.nmax=2, # use all if not specified
                     type = "LLO", # type of cv - LLO for now, after LTO, LLTO - CAST
                     k = 5, # number of folds
                     seed = 42,
                     folds, # if user want to create folds or column name
                     refit = TRUE,
                     output.format = "STFDF", # data.frame | sf | sftime | SpatVector      dodaj stars
                     parallel.processing = FALSE, # doParallel
                     pp.type = "snowfall", # "doParallel"
                     cpus=detectCores()-1,
                     progress=TRUE,
                     ...){
  
  # check the input
  message('Preparing data ...')
  if (missing(data) | missing(reg.coef) | missing(vgm.model)) {
    stop('The arguments data, reg.coef and vgm.model must not be empty!')
  }
  
  # prepare data
  if (class(data) %in% c("STFDF", "STSDF", "STIDF")) {
    if (is.numeric(obs.col)) {
      obs.col.name <- names(data@data)[obs.col]
    } else {
      obs.col.name <- obs.col
      obs.col <- index(names(data@data))[names(data@data) == obs.col.name]
    }
    if (is.numeric(data.staid.x.y.z[1])) {
      staid.name <- names(data@sp)[data.staid.x.y.z[1]]
    } else {
      staid.name <- data.staid.x.y.z[1]
      data.staid.x.y.z[1] <- index(names(data@sp))[names(data@sp) == staid.name]
    }
    s.crs <- data@sp@proj4string
  } else {
    if (!is.numeric(obs.col)) {
      obs.col.name <- obs.col
      obs.col <- index(names(data))[names(data) == obs.col.name]
    }
    data.prep <- data.prepare(data=data, data.staid.x.y.z=data.staid.x.y.z, obs.col=obs.col, s.crs=crs)
    data.df <- data.prep[["data.df"]]
    data.staid.x.y.z <- data.prep[["data.staid.x.y.z"]]
    staid.name <- names(data.df)[data.staid.x.y.z[1]]
    s.crs <- data.prep[["s.crs"]]
    if (is.na(s.crs)) {s.crs <- CRS(as.character(NA))}
    # obs.col.name <- data.prep[["obs.col"]]
    # to stfdf
    obs <- cbind(data.df[, c(data.staid.x.y.z[c(1,4)], obs.col)], data.df[, -c(data.staid.x.y.z[c(1:4)], obs.col)])
    if ('endTime' %in% names(obs)) { obs <- obs[, -which('endTime' == names(obs))] }
    if ('timeIndex' %in% names(obs)) { obs <- obs[, -which('timeIndex' == names(obs))] }
    stations <- data.df[, c(data.staid.x.y.z[c(1:3)])]
    stations <- unique(stations[complete.cases(stations), ])
    data <- meteo2STFDF(obs      = obs,
                        stations = stations,
                        crs      = s.crs, # CRS("+proj=longlat +datum=WGS84"),
                        obs.staid.time = c(1,2),
                        stations.staid.lon.lat = c(1,2,3)
    )
    obs.col=1
  }
  
  if (!inherits(data, "STFDF")) {
    data <- as(data, "STFDF")
  }
  
  # remove duplicates
  data <- rm.dupl(data, obs.col, zero.tol)
  
  names_covar <- names(reg.coef)[-1]
  # remove the stations where covariate is missing
  nrowsp <- length(data@sp)
  for (covar in names_covar){
    # count NAs per stations
    if (covar %in% names(data@data)) {
      numNA <- apply(matrix(data@data[,covar],
                            nrow=nrowsp,byrow= FALSE), MARGIN=1,
                     FUN=function(x) sum(is.na(x)))
      # Remove stations out of covariates
      rem <- numNA != length(time)
      data <-  data[rem,drop= FALSE]
    } else {
      data@sp <- data@sp[!is.na(data@sp@data[, covar]), ]
    }
  }
  
  time <- index(data@time)
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
  
  time <- index(data@time)
  data.df <- as.data.frame(data)
  
  c.dif <- setdiff(names_covar, names(data.df))
  
  # if covariate names doesn't exist in data - do ordinary kriging
  if (!identical(c.dif, character(0))) {
    # if covariates - do overlay
    warning(paste('The covariate(s) ', paste(c.dif, collapse = ", "), ' - missing from data!', sep = ""))
    warning('Trend is set to 0, performing space-time ordinary kriging cross-validation.')
    # data$tlm<-0 
    # data$tres <- data@data[,obs.col]
    names_covar <- names(data@data[obs.col])[1]
  }
  
  # set folds
  if (missing(folds)) {
    # create folds
    
    # check if time?
    
    # space_id <- rep(1:length(data@sp), length(time))
    # time_id <- rep(1:length(time), each = length(data@sp))
    # st_df <- cbind(space_id, time_id)
    # if (type == "LLO") {
    spacevar <- staid.name # names(data.df)[data.staid.x.y.z[1]] # "space_id"
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
    data$folds <- folds
    data.df$folds <- folds
    fold.column <- "folds"
  } else if (class(folds) %in% c("numeric", "character", "integer")) {
    # outer folds
    if (length(folds) == 1) { # column
      fold.column <- folds
      if (class(fold.column) %in% c("numeric", "integer")) {
        fold.column <- names(data)[fold.column]
      } else if (inherits(fold.column, "character")) {
        if (!fold.column %in% names(data)){
          stop(paste0('Colum with name "', fold.column, '" does not exist in data'))
        }
      }
      message(paste0("Fold column: ", fold.column))
    } else if (length(folds) == nrow(data.df)){ # vector
      data.df$folds <- folds
      fold.column <- "folds"
    } else {
      stop('length(folds) != nrow(data).')
    }
  }
  
  if(nrow(data.df[complete.cases(data.df), ]) == 0){
    stop('The data does not have complete cases!')
  } else {
    message(paste(nrow(data.df[complete.cases(data.df), ]), " observations complete cases (", length(data@sp), " stations x ", length(data@time), " days) used for cross-validation.", sep=""))
  }
  
  # doing CV
  message('Doing CV ...')
  pred <- c()
  vgm.model.init <- vgm.model
  for (val_fold in sort(unique(data.df[, fold.column]))) {
    message(paste('Fold ', val_fold, sep=""))
    
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
    
    fcv <- paste(obs.col.name, " ~ ", paste(names_covar, collapse = " + "), sep="")
    dev_df <- dev_df[, c(obs.col.name, names_covar)]
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
      dev$lm_res = dev@data[, obs.col] - dev$lm_trend
      
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
                           obs.col = obs.col.name,
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
    
    val_df <- as.data.frame(val)[, c(names(val_df)[1:4], staid.name, fold.column, obs.col.name)]
    val_df$pred <- val_stfdf$pred
    pred <- rbind(pred, val_df)
  }
  names(pred)[names(pred) == obs.col.name] <- "obs"
  
  # return
  if (output.format %in% c("sf", "sftime","SpatVector")) {
    sf = st_as_sf(pred, coords = c(1,2), crs = s.crs, agr = "constant") # sf
    if (output.format == "sf") {
      if (progress) message("Done!")
      return(sf)
    } else if (output.format == "sftime") {
      sftime <- st_sftime(sf)
      if (progress) message("Done!")
      return(sftime)
    } else if (output.format == "SpatVector") {
      sv <- vect(as(sf, "Spatial"))
      if (progress) message("Done!")
      return(sv)
    }
  } else if (output.format %in% c("STFDF", "STSDF", "STIDF")) {
    staid.index <- index(names(pred))[names(pred) == staid.name]
    sta <- pred[, c(staid.index,1,2)]
    sta <- sta[!duplicated(sta), ]
    obs <- pred[, c(staid.index, (4:length(pred))[which(4:length(pred) != staid.index)])] # 
    stfdf <- meteo2STFDF(obs      = obs,
                         stations = sta,
                         crs      = s.crs, # CRS("+proj=longlat +datum=WGS84"),
                         obs.staid.time = c(1,2),
                         stations.staid.lon.lat = c(1,2,3)
    )
    if (output.format == "STSDF") {
      if (progress) message("Done!")
      return(as(stfdf, "STSDF"))
    } else if (output.format == "STIDF") {
      if (progress) message("Done!")
      return(as(stfdf, "STIDF"))
    } else { #  (output.format == "STFDF")
      if (progress) message("Done!")
      return(stfdf)
    }
  } else {
    if (progress) message("Done!")
    return(pred)
  }
  ##########
}
