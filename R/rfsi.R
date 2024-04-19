rfsi <- function (formula, # without nearest obs
                  data, # sf | sftime | SpatVector | data.frame(x,y,obs,time,ec1,ec2,...)
                  data.staid.x.y.z = NULL, # = c(1,2,3,4), # if data.frame # time = mid.depth
                  n.obs = 5, # nearest obs
                  # time.nmax, # use all if not specified
                  avg = FALSE,
                  increment = 10000, # avg(nearest point dist)
                  range = 50000, # bbox smaller(a, b) / 2
                  quadrant = FALSE,
                  use.idw = FALSE,
                  idw.p = 2,
                  s.crs = NA, # to ce da leti - ...
                  p.crs = NA, # to ce da leti - ...
                  cpus = detectCores()-1, # for near.obs
                  progress = TRUE,
                  soil3d = FALSE, # soil RFSI
                  depth.range = 0.1, # in units of depth
                  no.obs = 'increase', # exactly
                  ...){ # ranger parameters
                        # quantreg,
                        # num.trees,
                        # mtry,
                        # min.node.size,
                        # sample.fraction,
  
  # check the input
  if (progress) message('Preparing data ...')
  if ((missing(data)) | missing(formula)) {
    stop('The arguments data (or obs and stations) and formula must not be empty!')
  }
  formula <- as.formula(formula)
  all_vars <- all.vars(formula)
  obs.col.name <- all_vars[1]
  names_covar <- all_vars[-1]
  
  # prepare data
  data.prep <- data.prepare(data=data, data.staid.x.y.z=data.staid.x.y.z, s.crs=s.crs)
  data.df <- data.prep[["data.df"]]
  data.staid.x.y.z <- data.prep[["data.staid.x.y.z"]]
  s.crs <- data.prep[["s.crs"]]
  
  # if vars exists
  c.dif <- setdiff(all_vars, names(data.df))
  if (!identical(c.dif, character(0))) {
    stop(paste('The variable(s) ', paste(c.dif, collapse = ", "), ' - missing from data!', sep = ""))
  }
  
  if (is.na(s.crs)) {
    warning('Source CRS is NULL! Using given coordinates for Euclidean distances calculation.')
  } else if (is.na(p.crs)) {
    if (progress) warning('Projection CRS is NULL. Using source CRS for Euclidean distances calculation:')
    if (progress) message(s.crs$input)
  } else if (identical(s.crs, p.crs)) {
    if (progress) message('s.crs==p.crs')
  } else {
    if (progress) message('Using projection CRS for Euclidean distances calculation.')
    prj_print <- try(paste('Do reprojection from: ', s.crs$input, ' to ', p.crs$input, sep=""), silent = TRUE)
    if(inherits(prj_print, "try-error")) {
      prj_print <- paste('Do reprojection from: ', s.crs@projargs, ' to ', p.crs@projargs, sep="")
    }
    if (progress) message(prj_print)
    # if (progress) message(paste('Do reprojection from: ', s.crs$input, ' to ', p.crs$input, sep=""))
    # reproject data coordinates
    data.coord <- data.df[, data.staid.x.y.z[2:3]]
    data.coord <- st_as_sf(data.coord, coords = names(data.coord), crs = s.crs, agr = "constant")
    data.coord <- st_transform(data.coord, p.crs)
    data.coord <- as.data.frame(st_coordinates(data.coord$geometry))
    names(data.coord) <- c('x.proj', 'y.proj')
    data.df <- cbind(data.df, data.coord)
    data.staid.x.y.z[2:3] <- c(length(data.df)-1, length(data.df))
  }
  
  x.y <- names(data.df)[data.staid.x.y.z[2:3]]
  data.df = data.df[complete.cases(data.df[, names_covar]), ]
  if (soil3d) {
    depth.name <- names(data.df)[data.staid.x.y.z[4]]
  }
  
  # if space-time
  if (!is.na(data.staid.x.y.z[4]) & !soil3d) {
    # sort data.df
    data.df <- data.df[order(data.df[, data.staid.x.y.z[4]],
                             data.df[, data.staid.x.y.z[1]]), ]
    time=sort(unique(data.df[, data.staid.x.y.z[4]]))
    daysNum = length(time)
    
    # calculate obs and dist
    if (progress) message('Space-time process ...')
    if (progress) message('Calculating distances to the nearest observations ...')
    registerDoParallel(cores=cpus)
    nearest_obs <- foreach (t = time, .export = c("near.obs")) %dopar% {
      
      day_df <- data.df[data.df[, data.staid.x.y.z[4]]==t, c(x.y, obs.col.name)]
      if (nrow(day_df)==0) {
        return(NULL)
      }
      return(near.obs(
        locations = day_df,
        # locations.x.y = c(1,2),
        observations = day_df,
        # observations.x.y = c(1,2),
        obs.col = obs.col.name,
        n.obs = n.obs,
        avg = avg,
        range = range,
        increment = increment,
        quadrant = quadrant,
        idw=use.idw,
        idw.p=idw.p
      ))
      
    }
    stopImplicitCluster()
    nearest_obs <- do.call("rbind", nearest_obs)
    
    # calculate TPS
    # if (use.tps) {
    #   if (progress) message('Calculating TPS ...')
    #   registerDoParallel(cores=cpus)
    #   tps_fit <- foreach (t = time) %dopar% {
    #     day_df <- data.df[data.df[, data.staid.x.y.z[4]]==t, c(x.y, obs.col.name)]
    #     if (nrow(day_df)==0) {
    #       return(NULL)
    #     }
    #     if ((nrow(day_df)-3) <= tps.df) {
    #       tps.df.day <- nrow(day_df)-3
    #     } else {
    #       tps.df.day <- tps.df
    #     }
    #     m <- Tps(day_df[, x.y], day_df[, obs.col.name],# lon.lat = T,
    #              #lambda = tps.lambda) #, GCV=F)
    #              df=tps.df.day)
    #     return(as.vector(m$fitted.values))
    #   }
    #   stopImplicitCluster()
    #   tps_fit <- unlist(tps_fit) #as.vector(do.call("rbind", tps_fit))
    # }
    
  } else if (!is.na(data.staid.x.y.z[4]) & soil3d) { # soil 3D
    # calculate obs and dist
    day_df <- data.df[, c(x.y, depth.name, obs.col.name)]
    if (progress) message('Soil 3D process ...')
    if (progress) message('Calculating distances to the nearest observations (profiles) ...')
    nearest_obs <- near.obs.soil(
      locations = day_df,
      # locations.x.y.md = c(1,2,3),
      observations = day_df,
      # observations.x.y.md = c(1,2,3),
      obs.col = 4,
      n.obs = n.obs,
      depth.range = depth.range,
      no.obs = no.obs,
      parallel.processing = TRUE,
      pp.type = "doParallel", # "snowfall"
      cpus = cpus
    )
    
  } else { # if spatial
    # calculate obs and dist
    day_df <- data.df[, c(x.y, obs.col.name)]
    if (progress) message('Spatial process ...')
    if (progress) message('Calculating distances to the nearest observations ...')
    nearest_obs <- near.obs(
      locations = day_df,
      # locations.x.y = c(1,2),
      observations = day_df,
      # observations.x.y = c(1,2),
      obs.col = obs.col.name,
      n.obs = n.obs,
      avg = avg,
      range = range,
      increment = increment,
      quadrant = quadrant,
      idw=use.idw,
      idw.p=idw.p
    )
    
    # # calculate TPS
    # if (use.tps) {
    #   if (progress) message('Calculating TPS ...')
    #   m <- Tps(day_df[, x.y], day_df[, obs.col.name],# lon.lat = T,
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

  data.df = data.df[complete.cases(data.df[, c(names_covar, names(nearest_obs))]), ]
  if (nrow(data.df) == 0) {
    stop("There is no complete cases in data.")
  }
  
  # if (use.tps) {
  #   formula = as.formula(paste(deparse(formula), paste(names(nearest_obs), collapse = " + "), tps.var, sep = " + "))
  # } else {
  formula = as.formula(paste(paste((deparse(formula)), collapse=""), paste(names(nearest_obs), collapse = " + "), sep = " + "))
  # }
  rm(nearest_obs)
  
  # fit RF model
  if (progress) message('Fitting RFSI model ...')
  rfsi_model <- ranger(formula, data = data.df, ...)
  # rfsi_model <- ranger(formula, data = data.df, importance = importance, seed = seed,
  # num.trees = num.trees, mtry = mtry, splitrule = "variance",
  # min.node.size = min.node.size,  sample.fraction = sample.fraction,
  # quantreg = quantreg)
  if (progress) message('Done!')
  
  return(rfsi_model)
  # dodaj formulu
  # dodaj n.obs
  # avg ...
  
}