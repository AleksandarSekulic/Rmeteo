pred.strk <- function (data, # data.frame(id,x,y,time,obs,ec1,ec2,...) | STFDF - with covariates
                       obs.col=1,
                       data.staid.x.y.z = NULL, # c(1,2,3,4), # if data.frame
                       newdata, # data.frame(x,y,time,ec1,ec2,...) | STFDF - with covariates
                       newdata.staid.x.y.z = NULL, # c(1,2,3), # if data.frame
                       z.value = NULL,
                       crs = NA, # brisi
                       zero.tol=0,
                       reg.coef, # check coef names
                       vgm.model,
                       sp.nmax=20, # use all if not specified
                       time.nmax=2, # use all if not specified
                       by='time', # 'station'
                       tiling= FALSE,
                       ntiles=64,
                       output.format = "STFDF", # data.frame | sf | sftime | SpatVector | SpatRaster      dodaj stars
                       parallel.processing = FALSE, # doParallel
                       pp.type = "snowfall", # "doParallel"
                       cpus=detectCores()-1,
                       computeVar=FALSE,
                       progress=TRUE,
                       ...){
  i <- NULL # to supress Warning
  # check the input
  if (progress) message('Preparing data ...')
  if ((missing(data)) | missing(newdata) | missing(reg.coef) | missing(vgm.model)) {
    stop('The arguments data, newdata, reg.coef and vgm.model must not be empty!')
  }
  
  # prepare data
  if (class(data) %in% c("STFDF", "STSDF", "STIDF")) {
    if (is.numeric(obs.col)) {
      obs.col.name <- names(data@data)[obs.col]
    } else {
      obs.col.name <- obs.col
      obs.col <- index(names(data@data))[names(data@data) == obs.col.name]
    }
  } else {
    if (!is.numeric(obs.col)) {
      obs.col.name <- obs.col
      obs.col <- index(names(data))[names(data) == obs.col.name]
    }
    data.prep <- data.prepare(data=data, data.staid.x.y.z=data.staid.x.y.z, obs.col=obs.col, s.crs=crs)
    data.df <- data.prep[["data.df"]]
    data.staid.x.y.z <- data.prep[["data.staid.x.y.z"]]
    s.crs <- data.prep[["s.crs"]]
    if (is.na(s.crs)) {s.crs <- CRS(as.character(NA))}
    obs.col.name <- data.prep[["obs.col"]]
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
  
  # prepare newdata
  if (class(newdata) %in% c("STFDF", "STSDF", "STIDF")) {
    newdata.s.crs <- newdata@sp@proj4string
  } else {
    newdata.prep <- data.prepare(data=newdata, data.staid.x.y.z=newdata.staid.x.y.z, s.crs=crs)
    newdata.df <- newdata.prep[["data.df"]]
    newdata.staid.x.y.z <- newdata.prep[["data.staid.x.y.z"]]
    newdata.s.crs <- newdata.prep[["s.crs"]]
    if (is.na(newdata.s.crs)) {newdata.s.crs <- CRS(as.character(NA))}
    # to stsdf
    obs <- cbind(newdata.df[, c(newdata.staid.x.y.z[c(1,4)])], newdata.df[, -c(newdata.staid.x.y.z[c(1:4)])])
    if ('endTime' %in% names(obs)) { obs <- obs[, -which('endTime' == names(obs))] }
    if ('timeIndex' %in% names(obs)) { obs <- obs[, -which('timeIndex' == names(obs))] }
    stations <- newdata.df[, c(newdata.staid.x.y.z[1:3])]
    stations <- unique(stations[complete.cases(stations), ])
    newdata <- meteo2STFDF(obs      = obs,
                           stations = stations,
                           crs      = newdata.s.crs, # CRS("+proj=longlat +datum=WGS84"),
                           obs.staid.time = c(1,2),
                           stations.staid.lon.lat = c(1,2,3)
    )
  }
  
  if (!inherits(data, "STFDF")) {
    data <- as(data, "STFDF")
  }
  if (!inherits(newdata, "STFDF")) {
    newdata <- as(newdata, "STFDF")
  }
  
  # check if data@time and newdata@time are of same class
  if (class(index(data@time)) != class(index(newdata@time))){
    stop(paste('The argument data and newdata must have time of the same class! ', class(index(data@time)), ' != ', class(index(newdata@time)) , sep=""))
  }
  
  names_covar <- names(reg.coef)[-1]
  
  # remove duplicates
  data <- rm.dupl(data, obs.col, zero.tol)
  newdata <- rm.dupl(newdata, 1, zero.tol)
  
  # DO OVERLAY !!!
  
  data.df <- as.data.frame(data)
  newdata.df <- as.data.frame(newdata)
  
  c.dif <- setdiff(names(reg.coef)[-1], names(newdata.df))
  if (!identical(c.dif, character(0))) {
    stop(paste('The covariate(s) ', paste(c.dif, collapse = ", "), ' - missing from newdata!', sep = ""))
  }
  newdata$tlm <- reg.coef[1] + as.matrix(newdata.df[, names(reg.coef)[-1]])  %*%  reg.coef[-1] # regression model-trend
  newdata$tlm <- as.vector(newdata$tlm)
  
  # remove the stations where covariate is missing
  nrowsp <- length(newdata@sp)
  for (covar in names_covar){
    # count NAs per stations
    if (covar %in% names(newdata@data)) {
      numNA <- apply(matrix(newdata@data[,covar],
                            nrow=nrowsp,byrow= FALSE), MARGIN=1,
                     FUN=function(x) sum(is.na(x)))
      # Remove stations out of covariates
      rem <- numNA != length(newdata@time)
      newdata <-  newdata[rem,drop= FALSE]
    } else {
      newdata@sp <- newdata@sp[!is.na(newdata@sp@data[, covar]), ]
    }
  }
  
  # Remove dates out of covariates
  rm.days <- c()
  for (t in 1:length(newdata@time)) {
    if(sum(complete.cases(newdata[, t]@data)) == 0) {
      rm.days <- c(rm.days, t)
    }
  }
  if(!is.null(rm.days)){
    newdata <- newdata[,-rm.days, drop= FALSE]
  }
  
  time <- newdata@time
  newdata.df <- as.data.frame(newdata)
  
  # newdata regression
  message('Doing regression ...')
  if(nrow(newdata.df[complete.cases(newdata.df), ]) == 0){
    warning('The argument newdata does not have complete cases! Trend is set to 0, performing space-time ordinary kriging.')
    newdata$tlm<-0 
    data$tres <- data@data[, obs.col]
  } else {
    # data regression
    c.dif <- setdiff(names(reg.coef)[-1], names(data.df))
    if (!identical(c.dif, character(0))) {
      warning(paste('The covariate(s) ', paste(c.dif, collapse = ", "), ' - missing from data! -> Doing overlay with newdata.', sep = ""))
      nrowsp <- length(data@sp)
      newdata@sp=as(newdata@sp,'SpatialPixelsDataFrame')
      ov <- sapply(1:length(time), function(i) over(data@sp, as(newdata[, i, 'tlm'], 'SpatialPixelsDataFrame')))
      ov <- do.call('cbind',ov)
      ov <- as.vector(ov)
      if (all(is.na(ov))) {
        stop(paste('There is no overlay of data with newdata!', sep = ""))
      }
      t1 <- which(as.character(index(time[1])) == as.character(index(data@time)))
      t2 <- which(as.character(index(time[length(time)])) == as.character(index(data@time)))
      data <- data[,t1:t2, drop= FALSE] # take only newdata days
      data$tlm <- ov
      
    } else {
      data$tlm <- reg.coef[1] + as.matrix(data.df[, names(reg.coef)[-1]])  %*%  reg.coef[-1] #regression model-trend
      data$tlm <- as.vector(data$tlm)
    }
    
    data$tres <- data@data[, obs.col]- data$tlm #residuals
    # count NAs per stations
    numNA <- apply(matrix(data@data[,'tres'],
                          nrow=nrowsp,byrow= FALSE), MARGIN=1,
                   FUN=function(x) sum(is.na(x)))
    # Remove stations out of covariates
    rem <- numNA != length(index(data@time))
    data <-  data[rem,drop= FALSE]
    
    # Remove dates out of covariates
    rm.days <- c()
    for (t in 1:length(index(data@time))) {
      if(sum(complete.cases(data[, t]@data)) == 0) {
        rm.days <- c(rm.days, t)
      }
    }
    if(!is.null(rm.days)){
      data <- data[,-rm.days, drop= FALSE]
    }
    
  } # end of regression
  
  # check sp.nmax
  nrowsp <- nrow(data@sp)
  if(sp.nmax > nrowsp) {
    warning(paste('The argument sp.nmax = ',sp.nmax , ' is larger than nrow(data@sp) = ', nrowsp, '! -> sp.nmax is set to ', nrowsp, '.', sep = ""))
    sp.nmax <- nrowsp
  }
  
  # # times - from to
  # i_1 <- (1:length(time)) - ceiling(time.nmax-1) # (index(time)) - ceiling(time.nmax/2)*24*60*60
  # i_1[i_1<1]=1        
  # ip1 <- i_1 + floor(time.nmax-1)
  # ip1[ip1>length(time)] <- length(time)
  
  # prediction
  message('Doing space-time kriging ...')
  message(paste('Doing for each loop by ', by, " ...", sep=""))
  newdata@sp=as(newdata@sp,'SpatialPointsDataFrame')
  newdata@sp$index=1:nrow(newdata@sp)
  # row.names(newdata@sp) = 1:nrow(newdata@sp)
  
  if (computeVar) {
    pred.var = 1:2
  } else {
    pred.var = 1
  }
  
  if(parallel.processing) {
    message(paste("Do parallel processing with", pp.type, "..."), sep="")
    if (pp.type == "doParallel") {
      registerDoParallel(cores=cpus)  
      cl <- makeCluster(cpus, type="SOCK")
    } else { # "snowfall"
      sfInit ( parallel = parallel.processing , cpus =cpus)
      sfLibrary(package="gstat", character.only=TRUE)
      sfLibrary(package="spacetime", character.only=TRUE)
      sfLibrary(package="sp", character.only=TRUE)
      sfExport("vgm.model" )
      sfExport( "pred.var" )
      # sfExport( "i_1" )
      # sfExport( "ip1" )
      sfExport( "time" )
      sfExport( "data" )
      sfExport( "newdata" )
      sfExport( "sp.nmax" )
    } 
  }
  
  if (!tiling) {
    temp.local<-data[,,'tres',drop= FALSE]
    if (progress)   pb <- txtProgressBar(style = 3,char= sprintf("pred krigeST ") , max=length(time) )
    if(parallel.processing & pp.type == "snowfall") {
      sfExport("temp.local" )
      sfExport( "progress" )
      sfExport( "pb" )
    }
    
    # for each by time or station
    if (by == 'station') {
      # function for kriging (because of doParallel)
      st=data@sp
      krige_fun <- function(i) {
        st$dist=spDists(temp.local@sp, newdata@sp[i, ])
        tmp_st<-st[ order(st$dist) ,]
        local_t= row.names(tmp_st[1:sp.nmax,])
        obs=temp.local[local_t, ,'tres', drop= FALSE]
        if (length(obs)<5){
          ret <- NA
        }
        # count NAs per stations
        numNA <- apply(matrix(obs@data[,'tres'],
                              nrow=length(obs@sp),byrow= FALSE), MARGIN=1,
                       FUN=function(x) sum(is.na(x)))
        # Remove stations out of covariates
        rem <- !numNA > 0
        obs <-  obs[rem,drop= FALSE]
        
        # If there are less than 5 observations
        if (length(obs)<5){
          ret <- NA
        } else {
          # nn <- newdata[i, , drop= FALSE]
          # index(nn@time) <- as.POSIXlt(index(nn@time))
          ret <- krigeST(as.formula("tres~1"),
                         data=as(obs, "STIDF"), 
                         newdata= newdata[1,, drop= FALSE],
                         modelList=vgm.model,
                         computeVar=T,
                         ...)@data[,pred.var]
          ret <- as.data.frame(ret)
          names(ret)[1] <- "var1.pred"
          ret$s_index <- i
          ret$t_index <- as.numeric(newdata@time)
        }
        if (progress)  setTxtProgressBar(pb, i )
        return(ret)
      }
      i_limit <- length(newdata@sp)
    } else {
      # function for kriging (because of doParallel)
      krige_fun <- function(i) {
        if (length(time) > 1) {
          sub_time <- index(newdata@time[i]) - (0:(time.nmax-1))*(index(newdata@time)[2] - index(newdata@time)[1])
        } else {
          sub_time <- index(newdata@time[i]) - (0:(time.nmax-1))*(newdata@endTime - as.POSIXlt(index(newdata@time)[1]))
        }
        ############# as.POSIXlt - maybe should be fixed !!!!!!!! #########
        
        sub_time <- sub_time[as.POSIXlt(sub_time) %in% as.POSIXlt(index(temp.local@time))]
        obs=temp.local[,as.POSIXlt(index(temp.local@time)) %in% as.POSIXlt(sub_time),'tres', drop= FALSE]
        if (length(obs)<5){
          ret <- NA
        }
        # count NAs per stations
        numNA <- apply(matrix(obs@data[,'tres'],
                              nrow=length(obs@sp),byrow= FALSE), MARGIN=1,
                       FUN=function(x) sum(is.na(x)))
        # Remove stations out of covariates
        rem <- !numNA > 0
        obs <-  obs[rem,drop= FALSE]

        # If there are less than 5 observations
        if (length(obs)<5){
          ret <- NA
        } else {
          ret <- krigeST(as.formula("tres~1"),
                         data=as(obs, "STIDF"),
                         newdata=STF(as(newdata@sp,"SpatialPoints"),
                                     temp.local@time[i],
                                     temp.local@endTime[i]),
                         modelList=vgm.model,
                         computeVar=T,
                         ...)@data[,pred.var]
          ret <- as.data.frame(ret)
          names(ret)[1] <- "var1.pred"
          ret$s_index <- newdata@sp$index
          ret$t_index <- i
        }
        if (progress)  setTxtProgressBar(pb, i )
        return(ret)
      }
      i_limit <- length(time)
    }
    
    if (parallel.processing) {
      if (pp.type == "doParallel") {
        xxx <- foreach(i = 1:i_limit, .packages = c("raster","spacetime","gstat")) %dopar% {krige_fun(i)}
      } else {
        xxx <- sfLapply (1:i_limit, function(i) {krige_fun(i)})
      }
    } else {
      xxx <- lapply (1:i_limit, function(i) {krige_fun(i)})
    }
    
    if (progress)  close(pb)
    res = do.call(rbind, xxx)
    res <- res[order(res$t_index, res$s_index), ]
    
  } else {
    # tiling
    message('Do tiling ...')
    dimnames(newdata@sp@coords)[[2]] <- c('x','y')
    xy=as.data.frame(newdata@sp)
    # xy= xy[row.names(newdata@sp),] # ???
    nxy= floor(sqrt(ntiles))
    xy$xg=as.character(cut(xy$x,nxy,labels=paste("x",1:nxy,sep="")))
    xy$yg=as.character(cut(xy$y,nxy,labels=paste("y",1:nxy,sep="")))
    xy$g=as.factor(paste(xy$xg,xy$yg,sep="") )
    xy$xg=NULL
    xy$yg=NULL
    coordinates(xy) = ~ x+y
    xy$index=1:nrow(xy)
    xy@proj4string <- newdata@sp@proj4string
    g_list <- split(xy, xy$g)
    # mAKE CHUNKS OF INITIAL DATA
    Mpoint=data.frame(x=mean(newdata@sp@coords[,1]),y=mean(newdata@sp@coords[,2]) )
    coordinates(Mpoint)=~x+y
    Mpoint@proj4string <- newdata@sp@proj4string
    st=data@sp
    # function for Middle point for each chunk (because of doParallel)
    mpts_fun <- function(i) {
      Mpoint=data.frame(x=mean(i@coords[,1]),y=mean(i@coords[,2]) )
      coordinates(Mpoint)=~x+y
      Mpoint@proj4string <- data@sp@proj4string
      Mpoint   }
    if (parallel.processing) {
      if (pp.type == "doParallel") {
        Mpts <- foreach(i = g_list, .packages = c("raster","spacetime","gstat")) %dopar% {mpts_fun(i)}
      } else {
        Mpts <- sfLapply (g_list, function(i) {mpts_fun(i)})
      }
    } else {
      Mpts <- lapply (g_list, function(i) {mpts_fun(i)})
    }
    temp.local <-data[ , ,'tres',drop= FALSE]
    
    if(parallel.processing & pp.type == "snowfall") {
      sfExport( "temp.local" )
      sfExport( "st" )
      sfExport( "Mpts" )
      sfExport( "g_list" ) }
    
    if (progress)   pb <- txtProgressBar(style = 3,char= sprintf("pred krigeST ") , max=length(g_list))
    if(parallel.processing & pp.type == "snowfall") {
      sfExport( "temp.local" )
      sfExport( "st" )
      sfExport( "Mpts" )
      sfExport( "g_list" )
      sfExport( "progress" )
      sfExport( "pb" )
    }
    # function for tile kriging (because of doParallel)
    krige_tile_fun <- function(i) {
      st$dist=spDists(temp.local@sp,Mpts[[i]])
      tmp_st<-st[ order(st$'dist') ,]
      local_t= row.names(tmp_st[1:sp.nmax,])
      
      xxx = as.list(rep(NA, length(time)))
      for( ii in 1:length(time) ) {
        sub_time <- index(newdata@time[i]) - (0:(time.nmax-1))*(index(newdata@time)[2] - index(newdata@time)[1])
        sub_time <- sub_time[as.POSIXlt(sub_time) %in% as.POSIXlt(index(temp.local@time))]
        # # obs = temp.local[local_t, i_1[ii]:ip1[ii],'tres',drop= FALSE]
        # obs = temp.local[, i_1[ii]:ip1[ii],'tres',drop= FALSE]
        obs=temp.local[local_t, as.POSIXlt(index(temp.local@time)) == as.POSIXlt(sub_time),'tres', drop= FALSE]
        if (length(obs)<5){
          xxx[[ii]] <- NA
        }

        # count NAs per stations
        numNA <- apply(matrix(obs@data[,'tres'],
                              nrow=length(obs@sp),byrow= FALSE), MARGIN=1,
                       FUN=function(x) sum(is.na(x)))
        # Remove stations out of covariates
        rem <- !numNA > 0
        obs <-  obs[rem,drop= FALSE]
        # If there are less than 5 observations
        if (length(obs)<5){
          xxx[[ii]] <- NA
        } else {
          xxx[[ii]] <- krigeST(as.formula("tres~1"),
                               data=as(obs, "STIDF"), 
                               newdata=STF(as(g_list[[i]],"SpatialPoints"),
                                           temp.local@time[ii],  
                                           temp.local@endTime[ii]),     
                               modelList=vgm.model,
                               # nmax = sp.nmax, ????????????????????????????
                               computeVar=T,
                               ...)@data[,pred.var]
          xxx[[ii]] <- as.data.frame(xxx[[ii]])
          names(xxx[[ii]])[1] <- "var1.pred"
          xxx[[ii]]$s_index <- g_list[[i]]$index
          xxx[[ii]]$t_index <- ii
        }
      } # end of  for
      if (progress)  setTxtProgressBar(pb, i )
      ret = do.call(rbind, xxx)
      ret }
    if (parallel.processing) {
      if (pp.type == "doParallel") {
        res <- foreach(i = 1:length(g_list), .packages = c("raster","spacetime","gstat")) %dopar% {krige_tile_fun(i)}
      } else {
        res <- sfLapply (1:length(g_list), function(i) {krige_tile_fun(i)})
      }
    } else {
      res <- lapply (1:length(g_list), function(i) {krige_tile_fun(i)})
    }
    if (progress)  close(pb)
    res = do.call(rbind, res)
    res <- res[order(res$t_index, res$s_index), ]
  } # end of tiling else do tiling
  
  # add results to stfdf
  resid <- as.numeric(res[, 1])
  newdata$pred <- resid + newdata$tlm
  stfdf <- newdata[ , , c("pred", "tlm"), drop= FALSE]
  if (computeVar) {
    p.var <- as.numeric(res[, 2])
    stfdf$var <- p.var
  }
  
  if (parallel.processing){
    if (pp.type == "doParallel"){
      stopImplicitCluster()
    } else {
      sfStop() 
    }
  }
  
  # return
  if (output.format == "STFDF") {
    if (progress) message("Done!")
    return(stfdf)
  } else if (output.format == "STIDF") {
    if (progress) message("Done!")
    return(as(stfdf, "STIDF"))
  } else if (output.format == "STSDF") {
    if (progress) message("Done!")
    return(as(stfdf, "STSDF"))
  } else if (output.format == "sftime") {
    sftime = st_as_sftime(as(stfdf, "STIDF"))
    if (progress) message("Done!")
    return(sftime)
  } else if (output.format == "sf") {
    sftime <- st_as_sftime(as(stfdf, "STIDF"))
    sf <- sftime
    sf <- st_drop_time(sf)
    sf$time <- sftime$time
    if (progress) message("Done!")
    return(sf)
  } else if (output.format == "SpatVector") {
    sftime <- st_as_sftime(as(stfdf, "STIDF"))
    sf <- sftime
    sf <- st_drop_time(sf)
    sf$time <- sftime$time
    sv <- vect(as(sf, "Spatial"))
    if (progress) message("Done!")
    return(sv)
  } else { # data.frame or SpatRaster
    result = stfdf
    result@sp <- as(result@sp, "SpatialPoints")
    result <- as.data.frame(as(result, "STIDF"))
    result <- result[, c(3,1:2,4, 7:(length(result)))]
    if (output.format == "data.frame") {
      if (progress) message("Done!")
      return(result)
    } else if (output.format == "SpatRaster") {
      if (!is.na(newdata.s.crs)) {newdata.s.crs <- st_crs(newdata.s.crs)$wkt}
      df <- result
      names(df)[4] <- "time"
      unique_times <- unique(df$time)
      f <- function(x) terra::rast(df[df$time == x, c(2:3,5:length(df))], type="xyz", crs = newdata.s.crs)
      sr <- sapply(unique_times, f)
      # sr <- rast(sr) # raster stack
      names(sr) <- unique_times # paste("pred_", unique_times, sep="")
      if (progress) message("Done!")
      return(sr)
    } else { #  (output.format == "data.frame")
      if (progress) message("Done!")
      return(result)
    }
  }
  ##########
}
  
  
  

