rfsi <- function (formula, # without nearest obs
                  data, # data.frame(x,y,obs,time,ec1,ec2,...) | STFDF - with covariates | SpatialPointsDataFrame | SpatialPixelsDataFrame
                  data.staid.x.y.time = c(1,2,3,4), # if data.frame
                  obs, # data.frame(id,time,obs,cov)
                  obs.staid.time = c(1,2),
                  stations, # data.frame(id,x,y)
                  stations.staid.x.y = c(1,2,3),
                  zero.tol = 0,
                  n.obs = 10, # nearest obs
                  # time.nmax, # use all if not specified
                  avg = FALSE,
                  increment = 10000, # avg(nearest point dist)
                  range = 50000, # bbox smaller(a, b) / 2
                  direct = FALSE,
                  use.idw = FALSE,
                  idw.p = 2,
                  s.crs = NA,
                  t.crs = NA,
                  cpus = detectCores()-1, # for near.obs
                  progress = TRUE,
                  ...){ # ranger parameters + quantreg!!!
  # num.trees,
  # mtry,
  # min.node.size,
  # sample.fraction,
  
  # check the input
  if (progress) print('Preparing data ...')
  if ((missing(data) & missing(obs) & missing(stations)) | missing(formula)) {
    stop('The arguments data (or obs and stations) and formula must not be empty!')
  }
  formula <- as.formula(formula)
  all_vars <- all.vars(formula)
  zcol.name <- all_vars[1]
  names_covar <- all_vars[-1]
  
  if (!missing(data)){
    if (class(data) == "data.frame") {
      # if data.staid.x.y.time is character
      if (!is.numeric(data.staid.x.y.time)) {
        data.staid.x.y.time <- sapply(data.staid.x.y.time, function(i) index(names(data))[names(data) == i])
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
      obs.staid.time <- sapply(obs.staid.time, function(i) index(names(obs))[names(obs) == i])
    }
    # if stations.staid.x.y is character
    if (!is.numeric(stations.staid.x.y)) {
      stations.staid.x.y <- sapply(stations.staid.x.y, function(i) index(names(stations))[names(stations) == i])
    }
    # to stfdf
    data.df <- join(obs, stations, by=names(obs)[obs.staid.time[1]], match="first")
    data.staid.x.y.time <- c(obs.staid.time[1],
                             stations.staid.x.y[2] + length(obs),
                             stations.staid.x.y[3] + length(obs),
                             obs.staid.time[2])
  }
  
  # if vars exists
  c.dif <- setdiff(all_vars, names(data.df))
  if (!identical(c.dif, character(0))) {
    stop(paste('The variable(s) ', paste(c.dif, collapse = ", "), ' - missing from data!', sep = ""))
  }
  
  if (is.na(s.crs)) {
    warning('Source CRS is NA! Using given coordinates for Euclidean distances calculation.')
  } else if (is.na(t.crs)) {
    if (progress) print('Targed CRS is NA. Using source CRS for Euclidean distances calculation:')
    if (progress) print(s.crs)
  } else {
    if (progress) print('Using target CRS for Euclidean distances calculation. Do reprojection from:')
    if (progress) print(s.crs)
    if (progress) print('to:')
    if (progress) print(t.crs)
    data.coord <- data.df[, data.staid.x.y.time[2:3]]
    coordinates(data.coord) <- names(data.coord)
    data.coord@proj4string <- s.crs
    data.coord <- spTransform(data.coord, t.crs)
    data.coord <- as.data.frame(data.coord)
    names(data.coord) <- c('x.proj', 'y.proj')
    data.df <- cbind(data.df, data.coord)
    data.staid.x.y.time[2:3] <- c(length(data.df)-1, length(data.df))
  }
  
  x.y <- names(data.df)[data.staid.x.y.time[2:3]]
  data.df = data.df[complete.cases(data.df), ]
  
  # if space-time
  if (!is.na(data.staid.x.y.time[4])) {
    # sort data.df
    data.df <- data.df[order(data.df[, data.staid.x.y.time[4]],
                             data.df[, data.staid.x.y.time[1]]), ]
    time=sort(unique(data.df[, data.staid.x.y.time[4]]))
    daysNum = length(time)
    
    # calculate obs and dist
    if (progress) print('Calculating distances to the nearest observations ...')
    registerDoParallel(cores=cpus)
    nearest_obs <- foreach (t = time, .export = c("near.obs")) %dopar% {
      
      day_df <- data.df[data.df[, data.staid.x.y.time[4]]==t, c(x.y, zcol.name)]
      if (nrow(day_df)==0) {
        return(NULL)
      }
      return(near.obs(
        locations = day_df,
        # locations.x.y = c(1,2),
        observations = day_df,
        # observations.x.y = c(1,2),
        zcol = zcol.name,
        n.obs = n.obs,
        avg = avg,
        range = range,
        increment = increment,
        direct = direct,
        idw=use.idw,
        idw.p=idw.p
      ))
      
    }
    stopImplicitCluster()
    nearest_obs <- do.call("rbind", nearest_obs)
    
    # calculate TPS
    # if (use.tps) {
    #   if (progress) print('Calculating TPS ...')
    #   registerDoParallel(cores=cpus)
    #   tps_fit <- foreach (t = time) %dopar% {
    #     day_df <- data.df[data.df[, data.staid.x.y.time[4]]==t, c(x.y, zcol.name)]
    #     if (nrow(day_df)==0) {
    #       return(NULL)
    #     }
    #     if ((nrow(day_df)-3) <= tps.df) {
    #       tps.df.day <- nrow(day_df)-3
    #     } else {
    #       tps.df.day <- tps.df
    #     }
    #     m <- Tps(day_df[, x.y], day_df[, zcol.name],# lon.lat = T,
    #              #lambda = tps.lambda) #, GCV=F)
    #              df=tps.df.day)
    #     return(as.vector(m$fitted.values))
    #   }
    #   stopImplicitCluster()
    #   tps_fit <- unlist(tps_fit) #as.vector(do.call("rbind", tps_fit))
    # }
    
  } else { # if spatial
    # calculate obs and dist
      day_df <- data.df[, c(x.y, zcol.name)]
      if (progress) print('Calculating distances to the nearest observations ...')
      nearest_obs <- near.obs(
        locations = day_df,
        # locations.x.y = c(1,2),
        observations = day_df,
        # observations.x.y = c(1,2),
        zcol = zcol.name,
        n.obs = n.obs,
        avg = avg,
        range = range,
        increment = increment,
        direct = direct,
        idw=use.idw,
        idw.p=idw.p
      )
      
      # # calculate TPS
      # if (use.tps) {
      #   if (progress) print('Calculating TPS ...')
      #   m <- Tps(day_df[, x.y], day_df[, zcol.name],# lon.lat = T,
      #            # lambda = tps.lambda) #, GCV=F)
      #            df=tps.df)
      #   tps_fit <- as.vector(m$fitted.values)
      # }
  }
  
  data.df <- cbind(data.df, nearest_obs)
  # if (use.tps) {
  #   tps.var <- paste("tps_", tps.df, sep="")
  #   data.df[, tps.var] <- tps_fit
  # }

  data.df = data.df[complete.cases(data.df), ]
  if (nrow(data.df) == 0) {
    stop("There is no complete cases in data.")
  }
  
  # if (use.tps) {
  #   formula = as.formula(paste(deparse(formula), paste(names(nearest_obs), collapse = " + "), tps.var, sep = " + "))
  # } else {
  formula = as.formula(paste(deparse(formula), paste(names(nearest_obs), collapse = " + "), sep = " + "))
  # }
  
  # fit RF model
  if (progress) print('Fitting RFSI model ...')
  rfsi_model <- ranger(formula, data = data.df, ...)
  # rfsi_model <- ranger(formula, data = data.df, importance = importance, seed = seed,
  # num.trees = num.trees, mtry = mtry, splitrule = "variance",
  # min.node.size = min.node.size,  sample.fraction = sample.fraction,
  # quantreg = quantreg)
  if (progress) print('Done!')
  
  return(rfsi_model)
  
}