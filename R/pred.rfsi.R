pred.rfsi <- function (model, # RFSI model
                       data, # sf | sftime | SpatVector | data.frame(x,y,z,obs)
                       obs.col=1,
                       data.staid.x.y.z = NULL, # if data.frame
                       newdata, # sf | sftime | SpatVector | SpatRaster | data.frame(x,y,z,ec1,ec2,...)
                       newdata.staid.x.y.z = NULL, # if data.frame
                       z.value = NULL,
                       zero.tol=0,
                       # time.nmax, # use all if not specified
                       s.crs = NA,
                       newdata.s.crs = NA,
                       p.crs = NA,
                       output.format = "data.frame",
                       cpus=detectCores()-1,
                       progress = TRUE,
                       soil3d = FALSE, # soil3D RFSI
                       depth.range = 0.1, # in units of depth
                       no.obs = 'increase', # exactly
                       ...){ # ranger parameters + quantiles!!!
  
  # check the input
  if (progress) print('Preparing data ...')
  if ((missing(data)) | missing(model) | missing(newdata)) {
    stop('The arguments model, data and newdata must not be empty!')
  }
  
  all_vars <- model$forest$independent.variable.names
  names_covar <- all_vars["obs" != substr(all_vars, 1, 3) & ("dist" != substr(all_vars, 1, 4) | all_vars=="dist") & "avg" != substr(all_vars, 1, 3)  & "dir" != substr(all_vars, 1, 3)  & "idw" != substr(all_vars, 1, 3)] # & "tps" != substr(all_vars, 1, 3) 
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
    quadrant = TRUE
  } else {
    quadrant = FALSE
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
  
  # prepare data
  data.prep <- data.prepare(data=data, data.staid.x.y.z=data.staid.x.y.z, obs.col=obs.col, s.crs=s.crs)
  data.df <- data.prep[["data.df"]]
  data.staid.x.y.z <- data.prep[["data.staid.x.y.z"]]
  s.crs <- data.prep[["s.crs"]]
  obs.col.name <- data.prep[["obs.col"]]
  
  # prepare newdata
  ############################################### terra #######################################
  # near.obs layeri u terra objekat
  # proveri names
  # predikcija direktno kroz ranger
  
  
  newdata.prep <- data.prepare(data=newdata, data.staid.x.y.z=newdata.staid.x.y.z, s.crs=newdata.s.crs)
  newdata.df <- newdata.prep[["data.df"]]
  newdata.staid.x.y.z <- newdata.prep[["data.staid.x.y.z"]]
  newdata.s.crs <- newdata.prep[["s.crs"]]
  
  # check data and new data domain
  st.class <- c("sftime") #, "data.frame")
  sp.class <- c("sf", "SpatVector", "data.frame")
  if (soil3d) { # soil3D - data.frame
    if (progress) print('Soil 3D process ...')
    data.cl <- "soil3D"
    if (!(any(class(data) %in% c(sp.class, "data.frame")))) {
      stop('The argument data must be or sf, SpatVector, or data.frame!')
    }
    if (is.na(data.staid.x.y.z[4])) {
      stop('The z (i.e. depth) element of argument data.staid.x.y.z must not be NA!')
    }
  } else { # not soil3D
    if (is.na(data.staid.x.y.z[4])) {
      if (progress) print('Spatial process ...')
      data.cl <- "s"
      data.staid.x.y.z <- data.staid.x.y.z[1:3]
    } else { # space-time
      if (progress) print('Space-time process ...')
      data.cl <- "st"
    }
  }
  
  sp.class <- c(sp.class, "SpatRaster")
  if (soil3d) { # soil3D - data.frame
    if (progress) print('Soil 3D process ...')
    newdata.cl <- "soil3D"
    if (!(any(class(newdata) %in% c(sp.class, "data.frame")))) {
      stop('The argument newdata must be sf, SpatVector, SpatRaster, or data.frame!')
    }
    if (any(class(newdata) == "SpatRaster")) {
      if (is.null(z.value) ) {
        stop('The arguments z.value must not be empty!')
      }
    } else {
      if (is.na(newdata.staid.x.y.z[4])) {
        stop('The z (i.e. depth) element of argument newdata.staid.x.y.z must not be NULL/NA!')
      }
    }
  } else { # not soil3D
    if (is.null(z.value) & is.na(newdata.staid.x.y.z[4])) {
      if (progress) print('Spatial process ...')
      newdata.cl <- "s"
      newdata.staid.x.y.z <- newdata.staid.x.y.z[1:3]
    } else { # space-time
      if (progress) print('Space-time process ...')
      newdata.cl <- "st"
      if (!(any(class(newdata) %in% c(sp.class, "data.frame")))) {
        stop('The argument newdata must be sf, SpatVector, SpatRaster, or data.frame!')
      }
      if (any(class(newdata) == "SpatRaster")) {
        if (is.null(z.value) ) {
          stop('The arguments z.value must not be empty!')
        }
      } else {
        if (is.na(newdata.staid.x.y.z[4])) {
          stop('The z (i.e. depth) element of argument newdata.staid.x.y.z must not be NULL/NA!')
        }
      }
    }
  }
  
  if (data.cl != newdata.cl) {
    stop('Arguments data and newdata must cover the same domain: spatial, spatio-temporal or soil3D!')
  }
  # check format
  if (newdata.cl %in% c("s", "soil3D") & output.format %in% c("sftime")) {
    stop('Argument output.format cannot be sftime for interpolation in spatial or soil3D domain!')
  }
  if (newdata.cl %in% c("st") & !(output.format %in% c("data.frame", "sftime", "SpatRaster"))) {
    stop('Argument output.format must be data.frame, sftime, or SpatRaster for interpolation in spatio-temporal domain!')
  }
  
  # if vars exists
  if (any(class(newdata) == "SpatRaster")) {
    c.dif <- setdiff(names_covar, names(newdata))
  } else {
    c.dif <- setdiff(names_covar, names(newdata.df))
  }
  if (!identical(c.dif, character(0))) {
    stop(paste('The variable(s) ', paste(c.dif, collapse = ", "), ' - missing from newdata!', sep = ""))
  }
  
  # check CRS
  if (is.na(s.crs)) {
    warning('Data source CRS is NA! Using given coordinates for Euclidean distances calculation.')
  } else if (is.na(p.crs)) {
    if (progress) print('Data projection CRS is NA Using source CRS for Euclidean distances calculation:')
    if (progress) print(s.crs$input)
  } else if (s.crs==p.crs) {
      if (progress) print('s.crs==p.crs')
  } else {
    if (progress) print('Using data projection CRS for Euclidean distances calculation.')
    if (progress) print(paste('Do reprojection from: ', s.crs$input, ' to ', p.crs$input, sep=""))
    # reproject data coordinates
    data.coord <- data.df[, data.staid.x.y.z[2:3]]
    data.coord <- st_as_sf(data.coord, coords = names(data.coord), crs = s.crs, agr = "constant")
    data.coord <- st_transform(data.coord, p.crs)
    data.coord <- as.data.frame(st_coordinates(data.coord$geometry))
    names(data.coord) <- c('x.proj', 'y.proj')
    data.df <- cbind(data.df, data.coord)
    data.staid.x.y.z[2:3] <- c(length(data.df)-1, length(data.df))
  }
  
  if (is.na(newdata.s.crs)) {
    warning('Newdata source CRS is NA! Using given coordinates for Euclidean distances calculation.')
    final.crs <- NA
  } else if (is.na(p.crs)) {
    if (progress) print('Newdata projection CRS is NA. Using source CRS for Euclidean distances calculation:')
    if (progress) print(newdata.s.crs$input)
    final.crs <- newdata.s.crs
  } else if (newdata.s.crs==p.crs) {
    if (progress) print('newdata.s.crs==p.crs')
    final.crs <- newdata.s.crs
  } else {
    old.newdata.x.y <- newdata.staid.x.y.z[2:3]
    if (progress) print('Using newdata projection CRS for Euclidean distances calculation.')
    if (progress) print(paste('Do reprojection from: ', s.crs$input, ' to ', p.crs$input, sep=""))
    # reproject data coordinates
    newdata.coord <- newdata.df[, newdata.staid.x.y.z[2:3]]
    newdata.coord <- st_as_sf(newdata.coord, coords = names(newdata.coord), crs = newdata.s.crs, agr = "constant")
    newdata.coord <- st_transform(newdata.coord, p.crs)
    newdata.coord <- as.data.frame(st_coordinates(newdata.coord$geometry))
    names(newdata.coord) <- c('x.proj', 'y.proj')
    newdata.df <- cbind(newdata.df, newdata.coord)
    newdata.staid.x.y.z[2:3] <- c(length(newdata.df)-1, length(newdata.df))
    final.crs <- newdata.s.crs
  }
  
  x.y <- names(data.df)[data.staid.x.y.z[2:3]]
  data.df = data.df[complete.cases(data.df[, (if(is.na(data.staid.x.y.z[4])) data.staid.x.y.z[-4] else data.staid.x.y.z)]), ]
  newdata.x.y <- names(newdata.df)[newdata.staid.x.y.z[2:3]]
  if (!any(class(newdata) == "SpatRaster")) {
    newdata.df = newdata.df[complete.cases(newdata.df[, names_covar]), ]
  }
  if (data.cl %in% c("st", "soil3D")) {
    data.3d.name <- names(data.df)[data.staid.x.y.z[4]]
    newdata.3d.name <- names(newdata.df)[newdata.staid.x.y.z[4]]
  }
  
  if (any(class(newdata) == "SpatRaster")) { # TERRA
    if (data.cl == "soil3D") { # soil 3D
      dev_day_df <- data.df[, c(x.y, data.3d.name, obs.col.name)]
      # if (!is.null(z.value)){
      rast.vect <- c()
      for (z.v in z.value) { # if length(z.value) > 1
        if (progress) print(paste('Z: ', z.v, sep=""))
        if (progress) print('Calculating distances to the nearest observations ...')
        # day_df <- cbind(crds(newdata.df, na.rm=FALSE, df=T), z.v)
        day_df <- newdata.df[, c(newdata.x.y, newdata.3d.name)]
        nearest_obs <- near.obs.soil(
          locations = day_df,
          # locations.x.y.md = c(1,2,3),
          observations = dev_day_df,
          # observations.x.y.md = c(1,2,3),
          obs.col = obs.col.name,
          n.obs = n.obs,
          depth.range = depth.range,
          no.obs = no.obs,
          parallel.processing = TRUE,
          pp.type = "doParallel", # "snowfall"
          cpus = cpus
        )
        # nlyr(newdata.df)
        r <- newdata.df[[1]]
        nearest_obs_list <- list()
        for (n in 1:length(names(nearest_obs))) {
          values(r) <- nearest_obs[, names(nearest_obs)[n]]
          names(r) <- names(nearest_obs)[n]
          nearest_obs_list[[n]] <- r
        }
        newdata.df.z <- c(newdata.df, rast(nearest_obs_list))
        
        # prediction
        if (progress) print('Doing RFSI predictions ...')
        pf <- function(model, ...) {
          library(ranger)
          predict(model, ..., num.threads = 1)$predictions
        }
        pred <- terra::predict(
          object = newdata.df.z,
          model = model,
          na.rm = TRUE,
          # type = "response" # "quantiles"
          # seed = 2021, # to control randomness
          # num.threads = 1, # to not overload RAM
          fun = pf,# function(model, ...) predict(model, ..., num.threads = 1)$predictions,
          cores = cpus,
          progress = "text",
          ...)
        if (progress) print('Newdata is in SpatRaster format. output.format ignored, returning SpatRaster!')
        if (length(names(pred)) == 1) {
          names(pred) <- "pred"
        } # else { # else quantiles
        #   for (n in 1:length(names(pred))) {
        #     names(pred)[n] <- paste(z.v, "_", names(pred)[n], sep="")
        #   }
        # }
        rast.vect <- c(rast.vect, pred)
      }
      names(rast.vect) <- z.value
      return(rast.vect)
      # } else { # reads from newdata.3d.name
      #   stop('The arguments z.value must not be empty!')
      #   day_df <- cbind(crds(newdata.df, na.rm=FALSE, df=T), values(newdata.df[[newdata.3d.name]]))
      #   nearest_obs <- near.obs.soil(
      #     locations = day_df,
      #     # locations.x.y.md = c(1,2,3),
      #     observations = dev_day_df,
      #     # observations.x.y.md = c(1,2,3),
      #     obs.col = obs.col.name,
      #     n.obs = n.obs,
      #     depth.range = depth.range,
      #     no.obs = no.obs,
      #     parallel.processing = TRUE,
      #     pp.type = "doParallel", # "snowfall"
      #     cpus = cpus
      #   )
      #   # nlyr(newdata.df)
      #   r <- newdata.df[[1]]
      #   nearest_obs_list <- list()
      #   for (n in 1:length(names(nearest_obs))) {
      #     values(r) <- nearest_obs[, names(nearest_obs)[n]]
      #     names(r) <- names(nearest_obs)[n]
      #     nearest_obs_list[[n]] <- r
      #   }
      #   newdata.df <- c(newdata.df, rast(nearest_obs_list))
      #   # prediction
      #   if (progress) print('Doing RFSI predictions ...')
      #   pf <- function(model, ...) {
      #     library(ranger)
      #     predict(model, ..., num.threads = 1)$predictions
      #   }
      #   pred <- terra::predict(
      #     object = newdata.df,
      #     model = model,
      #     na.rm = TRUE,
      #     # type = "response", # default
      #     # seed = 2021, # to control randomness
      #     # num.threads = 1, # to not overload RAM
      #     fun = pf,# function(model, ...) predict(model, ..., num.threads = 1)$predictions,
      #     cores = cpus,
      #     progress = "text",
      #     ...)
      #   if (progress) print('Newdata is in SpatRaster format. output.format ignored, returning SpatRaster!')
      #   if (length(names(pred)) == 1) {
      #     names(pred) <- paste("pred_", z.v, sep="")
      #   } else {# else quantiles
      #     for (n in 1:length(names(pred))) {
      #       names(pred)[n] <- paste(z.v, "_", names(pred)[n], sep="")
      #     }
      #   }
      #   return(pred)
      # }
    # end of soil 3D
    } else if (data.cl == "st") { # Space-time
      dev_day_df <- data.df[, c(x.y, data.3d.name, obs.col.name)]
      # if (!is.null(z.value)){
      rast.vect <- c()
      for (z.v in z.value) { # if length(z.value) > 1 ??? if there is no dynamic covariates
        if (progress) print(paste('Z: ', z.v, sep=""))
        if (progress) print('Calculating distances to the nearest observations ...')
        dev_day_df.z <- dev_day_df[dev_day_df[, 3] == z.v, ]
        # day_df <- cbind(crds(newdata.df, na.rm=FALSE, df=T), z.v)
        day_df <- newdata.df[, newdata.x.y]
        nearest_obs <- near.obs(
          locations = day_df,
          # locations.x.y = c(1,2),
          observations = dev_day_df.z,
          # observations.x.y = c(1,2),
          obs.col = obs.col.name,
          n.obs = n.obs,
          avg = avg,
          range = range,
          increment = increment,
          quadrant=quadrant,
          idw=use.idw,
          idw.p=idw.p
        )
        # nlyr(newdata.df)
        r <- newdata[[1]]
        nearest_obs_list <- list()
        for (n in 1:length(names(nearest_obs))) {
          values(r) <- nearest_obs[, names(nearest_obs)[n]]
          names(r) <- names(nearest_obs)[n]
          nearest_obs_list[[n]] <- r
        }
        newdata.z <- c(newdata, rast(nearest_obs_list))
        names(newdata.z)
        
        # prediction
        if (progress) print('Doing RFSI predictions ...')
        pf <- function(model, ...) {
          library(ranger)
          predict(model, ..., num.threads = 1)$predictions
        }
        pred <- terra::predict(
          object = newdata.z,
          model = model,
          na.rm = TRUE,
          # type = "response" # "quantiles"
          # seed = 2021, # to control randomness
          # num.threads = 1, # to not overload RAM
          fun = pf,# function(model, ...) predict(model, ..., num.threads = 1)$predictions,
          cores = cpus,
          progress = "text",
          ...)
        if (progress) print('Newdata is in SpatRaster format. output.format ignored, returning SpatRaster!')
        if (length(names(pred)) == 1) {
          names(pred) <- "pred"
        } # else { # else quantiles
        #   for (n in 1:length(names(pred))) {
        #     names(pred)[n] <- paste(z.v, "_", names(pred)[n], sep="")
        #   }
        # }
        rast.vect <- c(rast.vect, pred)
      }
      names(rast.vect) <- z.value
      return(rast.vect)
      # } else { # reads from newdata.staid.x.y.z[4]
      #   stop('The arguments z.value must not be empty!')
      # }
    } else { # if spatial
      # calculate obs and dist
      if (progress) print('Calculating distances to the nearest observations ...')
      dev_day_df <- data.df[, c(x.y, obs.col.name)]
      # day_df <- crds(newdata.df, na.rm=FALSE, df=T)
      day_df <- newdata.df[, newdata.x.y]
      nearest_obs <- near.obs(
        locations = day_df,
        # locations.x.y = c(1,2),
        observations = dev_day_df,
        # observations.x.y = c(1,2),
        obs.col = obs.col.name,
        n.obs = n.obs,
        avg = avg,
        range = range,
        increment = increment,
        quadrant=quadrant,
        idw=use.idw,
        idw.p=idw.p
      )
      # nlyr(newdata.df)
      r <- newdata[[1]]
      nearest_obs_list <- list()
      for (n in 1:length(names(nearest_obs))) {
        values(r) <- nearest_obs[, names(nearest_obs)[n]]
        names(r) <- names(nearest_obs)[n]
        nearest_obs_list[[n]] <- r
      }
      newdata <- c(newdata, rast(nearest_obs_list))
      # nlyr(newdata.df)
      
      # calculate TPS - TO DO
      # if (use.tps) {
      #   if (progress) print('Calculating TPS ...')
      #   m <- Tps(dev_day_df[, x.y], dev_day_df[, obs.col.name],# lon.lat = T,
      #            # lambda = tps.lambda) #, GCV=F)
      #            df=tps.df)
      #   tps_pred <- predict.Krig(m, day_df[, newdata.x.y])
      #   newdata.df[, tps.var] <- as.vector(tps_pred)
      # }
      
      # prediction
      if (progress) print('Doing RFSI predictions ...')
      pf <- function(model, ...) {
        library(ranger)
        predict(model, ..., num.threads = 1)$predictions
      }
      pred <- terra::predict(
        object = newdata,
        model = model,
        na.rm = TRUE,
        # type = "response" # "quantiles"
        # seed = 2021, # to control randomness
        # num.threads = 1, # to not overload RAM
        fun = pf,# function(model, ...) predict(model, ..., num.threads = 1)$predictions,
        cores = cpus,
        progress = "text",
        ...)
      if (progress) print('Newdata is in SpatRaster format. output.format ignored, returning SpatRaster!')
      if (length(names(pred)) == 1) {
        names(pred) <- "pred"
      } # else quantiles
      return(pred)
    }
    
  } else { # not SpatRaster
    # if space-time
    # if (!is.na(newdata.staid.x.y.z[4]) & !soil3d) {
    if (data.cl == "st") {
      time=sort(unique(newdata.df[, newdata.staid.x.y.z[4]]))
      daysNum = length(time)
      
      # calculate obs and dist
      if (progress) print('Calculating distances to the nearest observations ...')
      registerDoParallel(cores=cpus)
      nearest_obs <- foreach (t = time, .export = c("near.obs")) %dopar% {
        
        dev_day_df <- data.df[data.df[, data.staid.x.y.z[4]]==t, c(x.y, obs.col.name)]
        day_df <- newdata.df[newdata.df[, newdata.staid.x.y.z[4]]==t,
                             c(names(newdata.df)[newdata.staid.x.y.z[c(1,4)]],
                               newdata.x.y)]
        
        # day_df <- data.df[data.df$time==t, c(x.y, obs.col.name)]
        if (nrow(day_df)==0) {
          return(NULL)
        }
        ret <- cbind(day_df[, names(newdata.df)[newdata.staid.x.y.z[c(1,4)]]],
                     near.obs(
                       locations = day_df,
                       locations.x.y = c(3,4),
                       observations = dev_day_df,
                       # observations.x.y = c(1,2),
                       obs.col = obs.col.name,
                       n.obs = n.obs,
                       avg = avg,
                       range = range,
                       increment = increment,
                       quadrant=quadrant,
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
      #     dev_day_df <- data.df[data.df[, data.staid.x.y.z[4]]==t, c(x.y, obs.col.name)]
      #     day_df <- newdata.df[newdata.df[, newdata.staid.x.y.z[4]]==t,
      #                          c(names(newdata.df)[newdata.staid.x.y.z[c(1,4)]],
      #                            newdata.x.y)]
      #     if (nrow(day_df)==0) {
      #       return(NULL)
      #     }
      #     if ((nrow(dev_day_df)-3) <= tps.df) {
      #       tps.df.day <- nrow(dev_day_df)-3
      #     } else {
      #       tps.df.day <- tps.df
      #     }
      #     m <- Tps(dev_day_df[, x.y], dev_day_df[, obs.col.name],# lon.lat = T,
      #              # lambda = tps.lambda) #, GCV=F)
      #              df=tps.df.day)
      #     tps_pred <- cbind(day_df[, names(newdata.df)[newdata.staid.x.y.z[c(1,4)]]],
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
      
    } else if (data.cl == "soil3D") { # soil 3D
      if (progress) print('Calculating distances to the nearest observations ...')
      dev_day_df <- data.df[, c(x.y, data.3d.name, obs.col.name)]
      day_df <- newdata.df[, c(newdata.x.y, newdata.3d.name)]
      nearest_obs <- near.obs.soil(
        locations = day_df,
        # locations.x.y.md = c(1,2,3),
        observations = dev_day_df,
        # observations.x.y.md = c(1,2,3),
        obs.col = 4,
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
      dev_day_df <- data.df[, c(x.y, obs.col.name)]
      day_df <- newdata.df[, c(newdata.x.y)]
      nearest_obs <- near.obs(
        locations = day_df,
        # locations.x.y = c(1,2),
        observations = dev_day_df,
        # observations.x.y = c(1,2),
        obs.col = obs.col.name,
        n.obs = n.obs,
        avg = avg,
        range = range,
        increment = increment,
        quadrant=quadrant,
        idw=use.idw,
        idw.p=idw.p
      )
      newdata.df <- cbind(newdata.df, nearest_obs)
      
      # calculate TPS
      # if (use.tps) {
      #   if (progress) print('Calculating TPS ...')
      #   m <- Tps(dev_day_df[, x.y], dev_day_df[, obs.col.name],# lon.lat = T,
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
    
    # prediction
    if (progress) print('Doing RFSI predictions ...')
    pred <- predict(model, newdata.df, ...)$predictions
    if (exists("old.newdata.x.y")) {
      if (length(newdata.staid.x.y.z) == 4) {
        result <- cbind(newdata.df[, c(newdata.staid.x.y.z[1], old.newdata.x.y, newdata.staid.x.y.z[c(4,2,3)])], pred)
      } else {
        result <- cbind(newdata.df[, c(newdata.staid.x.y.z[1], old.newdata.x.y, newdata.staid.x.y.z[c(2,3)])], pred)
      }
    } else {
      result <- cbind(newdata.df[, newdata.staid.x.y.z], pred)
    }
  
  }
  
  # return
  if (output.format == "sftime") {
    sf <- result
    names(sf)[4] <- "time"
    sf <- st_as_sf(sf, coords = names(sf)[2:3], crs = final.crs, agr = "constant")
    # sf$time <- Sys.time()
    sf$time <- as.POSIXct(sf$time)
    sftime <- st_sftime(sf)
    if (progress) print("Done!")
    return(sftime)
  } else if (output.format == "sf") {
    sf <- st_as_sf(result, coords = names(result)[2:3], crs = final.crs, agr = "constant")
    if (progress) print("Done!")
    return(sf)
  } else if (output.format == "SpatVector") {
    if (!is.na(final.crs)) {final.crs <- final.crs$wkt}
    sv <- terra::vect(result, geom = names(result)[2:3], crs = final.crs)
    if (progress) print("Done!")
    return(sv)
  } else if (output.format == "SpatRaster") {
    if (!is.na(final.crs)) {final.crs <- final.crs$wkt}
    if (newdata.cl == "st"){
      df <- result
      names(df)[4] <- "time"
      unique_times <- unique(df$time)
      if (exists("old.newdata.x.y")) {
        f <- function(x) terra::rast(df[df$time == x, c(2:3,7:length(df))], type="xyz", crs = final.crs)
      } else {
        f <- function(x) terra::rast(df[df$time == x, c(2:3,5:length(df))], type="xyz", crs = final.crs)
      }
      sr <- sapply(unique_times, f)
      # sr <- rast(sr) # raster stack
      # names(sr) <- paste("pred_", unique_times, sep="")
    } else {
      if (exists("old.newdata.x.y")) {
        sr <- terra::rast(result[, c(2:3,6:length(result))], type="xyz", crs = final.crs)
      } else {
        sr <- terra::rast(result[, c(2:3,4:length(result))], type="xyz", crs = final.crs)
      }
    }
    if (progress) print("Done!")
    return(sr)
  } else { #  (output.format == "df")
    if (progress) print("Done!")
    return(result)
  }
  
}

# data.frame - + pred, pred.qXX, ...
# STFDF - + pred, pred.qXX, ...