pred.rfsi <- function (model, # RFSI model iz rfsi ili tune rfsi funkcije
                       data, # data.frame(x,y,obs,time) | STFDF - with covariates | SpatialPointsDataFrame | SpatialPixelsDataFrame
                       zcol=1,
                       data.staid.x.y.time = c(1,2,3,4), # if data.frame
                       obs, # data.frame(id,time,obs,cov)
                       obs.staid.time = c(1,2),
                       stations, # data.frame(id,x,y)
                       stations.staid.x.y = c(1,2,3),
                       newdata, # data.frame(x,y,time,ec1,ec2,...) | STFDF - with covariates | SpatialPointsDataFrame | SpatialPixelsDataFrame
                       newdata.staid.x.y.time = c(1,2,3), # if data.frame
                       zero.tol=0,
                       # time.nmax, # use all if not specified
                       s.crs=NA,
                       newdata.s.crs=NA,
                       t.crs=NA,
                       output.format = "data.frame",
                       cpus=detectCores()-1,
                       progress=TRUE,
                       soil3d = FALSE, # soil RFSI
                       depth.range = 0.1, # in units of depth
                       no.obs = 'increase', # exactly
                       ...){ # ranger parameters + quantreg!!! - ovoga nema!!!
  
  # check the input
  if (progress) print('Preparing data ...')
  if ((missing(data) & missing(obs) & missing(stations)) | missing(model) | missing(newdata)) {
    stop('The arguments model, data (or obs and stations) and newdata must not be empty!')
  }
  
  all_vars <- model$forest$independent.variable.names
  names_covar <- all_vars["obs" != substr(all_vars, 1, 3) & "dist" != substr(all_vars, 1, 4) & "avg" != substr(all_vars, 1, 3)  & "dir" != substr(all_vars, 1, 3)  & "idw" != substr(all_vars, 1, 3)] # & "tps" != substr(all_vars, 1, 3) 
  n.obs <- sum("obs" == substr(model$forest$independent.variable.names, 1, 3))
  if("avg" %in% substr(model$forest$independent.variable.names, 1, 3)){
    avgs <- model$forest$independent.variable.names["avg" == substr(model$forest$independent.variable.names, 1, 3)]
    avg = TRUE
    increment <- as.numeric(substr(avgs[1], 4, nchar(avgs[1])))
    range <- as.numeric(substr(avgs[length(avgs)], 4, nchar(avgs[length(avgs)])))
  } else {
    avg = FALSE
    increment <- NULL
    range <- NULL
  }
  
  if("dir" %in% substr(model$forest$independent.variable.names, 1, 3)){
    direct = TRUE
  } else {
    direct = FALSE
  }
  # if("tps" %in% substr(model$forest$independent.variable.names, 1, 3)){
  #   use.tps = TRUE
  #   tps.var <- model$forest$independent.variable.names["tps" == substr(model$forest$independent.variable.names, 1, 3)][1]
  #   tps.df = as.numeric(strsplit(tps.var, "_")[[1]][2])
  # } else {
  #   use.tps = FALSE
  # }
  
  if("idw" %in% substr(model$forest$independent.variable.names, 1, 3)){
    use.idw = TRUE
    idw.var <- model$forest$independent.variable.names["idw" == substr(model$forest$independent.variable.names, 1, 3)]#[1]
    # idw.p = as.numeric(strsplit(idw.var, "_")[[1]][2])
    idw.p = as.numeric(sapply(idw.var, function(x) strsplit(x, "_")[[1]][2]))
  } else {
    use.idw = FALSE
  }
  
  if (!missing(data)){
    if (class(data) == "data.frame") {
      # if zcol is character
      if (is.numeric(zcol)) {
        zcol.name <- names(data)[zcol]
      } else {
        zcol.name <- zcol
      }
      if (!is.numeric(data.staid.x.y.time)) {
        data.staid.x.y.time <- match(data.staid.x.y.time, names(data)) # sapply(data.staid.x.y.time, function(i) index(names(data))[names(data) == i])
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
      if (is.numeric(zcol)) {
        zcol.name <- names(data@data)[zcol]
      } else {
        zcol.name <- zcol
      }
    } else if (class(data) == "SpatialPointsDataFrame" | class(data) == "SpatialPixelsDataFrame") {
      data.df <- as.data.frame(data)
      data.df$staid <- 1:nrow(data.df)
      data.staid.x.y.time <- c(length(data.df),length(data.df)-2,length(data.df)-1,NA)
      if (!is.na(data@proj4string)) {
        s.crs <- data@proj4string
      }
      if (is.numeric(zcol)) {
        zcol.name <- names(data@data)[zcol]
      } else {
        zcol.name <- zcol
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
    data.df <- join(obs, stations, by=names(obs)[obs.staid.time[1]], match="first")
    data.staid.x.y.time <- c(obs.staid.time[1],
                             stations.staid.x.y[2] + length(obs),
                             stations.staid.x.y[3] + length(obs),
                             obs.staid.time[2])
    if (is.numeric(zcol)) {
      zcol.name <- names(obs)[zcol]
    } else {
      zcol.name <- zcol
    }
  }
  
  # newdata
  if (class(newdata) == "data.frame") {
    if (!is.numeric(newdata.staid.x.y.time)) {
      newdata.staid.x.y.time <- match(newdata.staid.x.y.time, names(newdata)) # sapply(data.staid.x.y.time, function(i) index(names(data))[names(data) == i])
    }
    newdata.df = newdata
  } else if (class(newdata) == "STFDF" | class(newdata) == "STSDF") {
    if (class(newdata) == "STSDF") {newdata <- as(newdata, "STFDF")}
    newdata <- rm.dupl(newdata, 1, zero.tol)
    newdata.df <- as.data.frame(newdata)
    newdata.staid.x.y.time <- c(3,1,2,4)
    if (!is.na(newdata@sp@proj4string)) {
      newdata.s.crs <- newdata@sp@proj4string
    }
  } else if (class(newdata) == "SpatialPointsDataFrame" | class(newdata) == "SpatialPixelsDataFrame") {
    newdata.df <- as.data.frame(newdata)
    newdata.df$staid <- 1:nrow(newdata.df)
    newdata.staid.x.y.time <- c(length(newdata.df),length(newdata.df)-2,length(newdata.df)-1,NA)
    if (!is.na(newdata@proj4string)) {
      newdata.s.crs <- newdata@proj4string
    }
  } else {
    stop('The argument newdata must be of STFDF, STSDF, data.frame,SpatialPointsDataFrame or SpatialPixelsDataFrame class!') # "STSDF"
  }
  
  # check data and new data domain
  st.class <- c("STFDF", "STSDF", "data.frame")
  if (class(data) %in% st.class & !is.na(data.staid.x.y.time[4]) & !soil3d) { # space-time
    if (progress) print('Space-time process ...')
    data.cl <- "st"
  } else if (!(class(data) %in% st.class[1:2]) & !is.na(data.staid.x.y.time[4]) & soil3d) { # soil
    if (progress) print('Soil 3D process ...')
    data.cl <- "soil"
  } else { # spatial
    if (progress) print('Spatial process ...')
    data.cl <- "s"
    data.staid.x.y.time <- data.staid.x.y.time[1:3]
  }
  if (class(newdata) %in% st.class & !is.na(newdata.staid.x.y.time[4]) & !soil3d) {
    newdata.cl <- "st"
  } else if (!(class(newdata) %in% st.class[1:2]) & !is.na(newdata.staid.x.y.time[4]) & soil3d) { # soil
    newdata.cl <- "soil"
  } else {
    newdata.cl <- "s"
    newdata.staid.x.y.time <- newdata.staid.x.y.time[1:3]
  }
  if (data.cl != newdata.cl) {
    stop('Arguments data and newdata must cover the same domain: spatial, spatio-temporal or soil!')
  }
  # check format
  if (newdata.cl %in% c("s", "soil") & output.format %in% c("STFDF", "STSDF")) {
    stop('Argument output.format cannot be STFDF or STSDF for interpolation in spatial or soil domain!')
  }
  
  # if vars exists
  c.dif <- setdiff(names_covar, names(newdata.df))
  if (!identical(c.dif, character(0))) {
    stop(paste('The variable(s) ', paste(c.dif, collapse = ", "), ' - missing from newdata!', sep = ""))
  }
  
  # check CRS
  if (is.na(s.crs)) {
    warning('Data source CRS is NA! Using given coordinates for Euclidean distances calculation.')
  } else if (is.na(t.crs)) {
    if (progress) print('Data targed CRS is NA. Using source CRS for Euclidean distances calculation:')
    if (progress) print(s.crs)
  } else {
    if (progress) print('Using data target CRS for Euclidean distances calculation. Do reprojection from:')
    if (progress) print(s.crs)
    if (progress) print('to:')
    if (progress) print(t.crs)
    # reproject data coordinates
    data.coord <- data.df[, data.staid.x.y.time[2:3]]
    coordinates(data.coord) <- names(data.coord)
    data.coord@proj4string <- s.crs
    data.coord <- spTransform(data.coord, t.crs)
    data.coord <- as.data.frame(data.coord)
    names(data.coord) <- c('x.proj', 'y.proj')
    data.df <- cbind(data.df, data.coord)
    data.staid.x.y.time[2:3] <- c(length(data.df)-1, length(data.df))
  }
  
  if (is.na(newdata.s.crs)) {
    warning('Newdata source CRS is NA! Using given coordinates for Euclidean distances calculation.')
    final.crs <- NA
  } else if (is.na(t.crs)) {
    if (progress) print('Newdata targed CRS is NA. Using source CRS for Euclidean distances calculation:')
    if (progress) print(newdata.s.crs)
    final.crs <- newdata.s.crs
  } else {
    old.newdata.x.y <- newdata.staid.x.y.time[2:3]
    if (progress) print('Using newdata target CRS for Euclidean distances calculation. Do reprojection from:')
    if (progress) print(newdata.s.crs)
    if (progress) print('to:')
    if (progress) print(t.crs)
    # reproject data coordinates
    newdata.coord <- newdata.df[, newdata.staid.x.y.time[2:3]]
    coordinates(newdata.coord) <- names(newdata.coord)
    newdata.coord@proj4string <- newdata.s.crs
    newdata.coord <- spTransform(newdata.coord, t.crs)
    newdata.coord <- as.data.frame(newdata.coord)
    names(newdata.coord) <- c('x.proj', 'y.proj')
    newdata.df <- cbind(newdata.df, newdata.coord)
    newdata.staid.x.y.time[2:3] <- c(length(newdata.df)-1, length(newdata.df))
    final.crs <- newdata.s.crs
  }
  
  x.y <- names(data.df)[data.staid.x.y.time[2:3]]
  data.df = data.df[complete.cases(data.df[, (if(is.na(data.staid.x.y.time[4])) data.staid.x.y.time[-4] else data.staid.x.y.time))]), ]
  newdata.x.y <- names(newdata.df)[newdata.staid.x.y.time[2:3]]
  newdata.df = newdata.df[complete.cases(newdata.df[, names_covar]), ]
  if (soil3d) {
    data.depth.name <- names(data.df)[data.staid.x.y.time[4]]
    newdata.depth.name <- names(newdata.df)[newdata.staid.x.y.time[4]]
  }
  
  # if space-time
  # if (!is.na(newdata.staid.x.y.time[4]) & !soil3d) {
  if (data.cl == "st") {
    time=sort(unique(newdata.df[, newdata.staid.x.y.time[4]]))
    daysNum = length(time)
    
    # calculate obs and dist
    if (progress) print('Calculating distances to the nearest observations ...')
    registerDoParallel(cores=cpus)
    nearest_obs <- foreach (t = time, .export = c("near.obs")) %dopar% {
      
      dev_day_df <- data.df[data.df[, data.staid.x.y.time[4]]==t, c(x.y, zcol.name)]
      day_df <- newdata.df[newdata.df[, newdata.staid.x.y.time[4]]==t,
                           c(names(newdata.df)[newdata.staid.x.y.time[c(1,4)]],
                             newdata.x.y)]
      
      # day_df <- data.df[data.df$time==t, c(x.y, zcol.name)]
      if (nrow(day_df)==0) {
        return(NULL)
      }
      ret <- cbind(day_df[, names(newdata.df)[newdata.staid.x.y.time[c(1,4)]]],
                   near.obs(
                     locations = day_df,
                     locations.x.y = c(3,4),
                     observations = dev_day_df,
                     # observations.x.y = c(1,2),
                     zcol = zcol.name,
                     n.obs = n.obs,
                     avg = avg,
                     range = range,
                     increment = increment,
                     direct=direct,
                     idw=use.idw,
                     idw.p=idw.p
                   ))
      return(ret)
    }
    stopImplicitCluster()
    nearest_obs <- do.call("rbind", nearest_obs)
    # join by staid, date
    newdata.df <- join(newdata.df, nearest_obs)# cbind(data.df, nearest_obs)
    
    # calculate TPS
    # if (use.tps) {
    #   if (progress) print('Calculating TPS ...')
    #   registerDoParallel(cores=cpus)
    #   tps_fit <- foreach (t = time) %dopar% {
    #     dev_day_df <- data.df[data.df[, data.staid.x.y.time[4]]==t, c(x.y, zcol.name)]
    #     day_df <- newdata.df[newdata.df[, newdata.staid.x.y.time[4]]==t,
    #                          c(names(newdata.df)[newdata.staid.x.y.time[c(1,4)]],
    #                            newdata.x.y)]
    #     if (nrow(day_df)==0) {
    #       return(NULL)
    #     }
    #     if ((nrow(dev_day_df)-3) <= tps.df) {
    #       tps.df.day <- nrow(dev_day_df)-3
    #     } else {
    #       tps.df.day <- tps.df
    #     }
    #     m <- Tps(dev_day_df[, x.y], dev_day_df[, zcol.name],# lon.lat = T,
    #              # lambda = tps.lambda) #, GCV=F)
    #              df=tps.df.day)
    #     tps_pred <- cbind(day_df[, names(newdata.df)[newdata.staid.x.y.time[c(1,4)]]],
    #                       predict.Krig(m, day_df[, newdata.x.y]))
    #     names(tps_pred)[3] <- tps.var
    #     return(tps_pred)
    #   }
    #   stopImplicitCluster()
    #   
    #   tps_fit <- do.call("rbind", tps_fit)
    #   # join by staid, date
    #   newdata.df <- join(newdata.df, tps_fit)# cbind(data.df, nearest_obs)
    # }
    
  } else if (data.cl == "soil") { # soil 3D
    if (progress) print('Calculating distances to the nearest observations ...')
    dev_day_df <- data.df[, c(x.y, data.depth.name, zcol.name)]
    day_df <- newdata.df[, c(newdata.x.y, newdata.depth.name)]
    nearest_obs <- near.obs.soil(
      locations = day_df,
      # locations.x.y.md = c(1,2,3),
      observations = dev_day_df,
      # observations.x.y.md = c(1,2,3),
      zcol = 4,
      n.obs = n.obs,
      depth.range = depth.range,
      no.obs = no.obs,
      parallel.processing = TRUE,
      pp.type = "doParallel", # "snowfall"
      cpus = cpus
    )
    newdata.df <- cbind(newdata.df, nearest_obs)
  } else { # if spatial
    # calculate obs and dist
    if (progress) print('Calculating distances to the nearest observations ...')
    dev_day_df <- data.df[, c(x.y, zcol.name)]
    day_df <- newdata.df[, c(newdata.x.y)]
    nearest_obs <- near.obs(
      locations = day_df,
      # locations.x.y = c(1,2),
      observations = dev_day_df,
      # observations.x.y = c(1,2),
      zcol = zcol.name,
      n.obs = n.obs,
      avg = avg,
      range = range,
      increment = increment,
      direct=direct,
      idw=use.idw,
      idw.p=idw.p
    )
    newdata.df <- cbind(newdata.df, nearest_obs)
    
    # calculate TPS
    # if (use.tps) {
    #   if (progress) print('Calculating TPS ...')
    #   m <- Tps(dev_day_df[, x.y], dev_day_df[, zcol.name],# lon.lat = T,
    #            # lambda = tps.lambda) #, GCV=F)
    #            df=tps.df)
    #   tps_pred <- predict.Krig(m, day_df[, newdata.x.y])
    #   newdata.df[, tps.var] <- as.vector(tps_pred)
    # }
    
  }
  
  newdata.df = newdata.df[complete.cases(newdata.df[, c(names_covar, names(nearest_obs))]), ]
  if (nrow(newdata.df) == 0) {
    stop("There is no complete cases in newdata.")
  }
  
  if (progress) print('Making RFSI predictions ...')
  pred <- predict(model, newdata.df, ...)$predictions
  if (exists("old.newdata.x.y")) {
    if (length(newdata.staid.x.y.time) == 4) {
      result <- cbind(newdata.df[, c(newdata.staid.x.y.time[1], old.newdata.x.y, newdata.staid.x.y.time[c(4,2,3)])], pred)
    } else {
      result <- cbind(newdata.df[, c(newdata.staid.x.y.time[1], old.newdata.x.y, newdata.staid.x.y.time[c(2,3)])], pred)
    }
  } else {
    result <- cbind(newdata.df[, newdata.staid.x.y.time], pred)
  }
  
  # return
  if (output.format == "STSDF" | output.format == "STFDF" ) {
    if (exists("old.newdata.x.y")) {
      stations <- result[, c(1,2:3,5:6)]
      obs <- result[, c(1,4,7:length(result))]
    } else {
      stations <- result[, 1:3]
      obs <- result[, c(1,4:length(result))]
    }
    stations <- unique(stations[complete.cases(stations), ])
    stfdf <- meteo2STFDF(obs      = obs,
                         stations = stations,
                         crs      = final.crs,
                         obs.staid.time = c(1,2),
                         stations.staid.lon.lat = c(1,2,3))
    if (output.format == "STSDF") {
      stfdf <- as(stfdf, "STSDF")
    }
    if (progress) print("Done!")
    return(stfdf)
  } else if (output.format == "SpatialPointsDataFrame") {
    spdf <- result
    coordinates(spdf) <- c(names(spdf)[2:3])
    if(!is.na(final.crs)){
      spdf@proj4string <- final.crs
    }
    if (progress) print("Done!")
    return(spdf)
  } else if (output.format == "SpatialPixelsDataFrame") {
    spdf <- result
    coordinates(spdf) <- c(names(spdf)[2:3])
    if(!is.na(final.crs)){
      spdf@proj4string <- final.crs
    }
    spdf <- as(spdf, "SpatialPixelsDataFrame")
    if (progress) print("Done!")
    return(spdf)
  } else { #  (output.format == "df")
    if (progress) print("Done!")
    return(result)
  }
  
}

# data.frame - + pred, pred.qXX, ...
# STFDF - + pred, pred.qXX, ...