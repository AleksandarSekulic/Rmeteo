meteo2STFDF <- function(obs,
                        stations,
                        obs.staid.time=c(1,2),
                        stations.staid.lon.lat=c(1,2,3),
                        crs=CRS(as.character(NA)),
                        delta=NULL
                        ) {
 
  if (class(stations.staid.lon.lat) == "character") {
    stations.staid.lon.lat <- match(stations.staid.lon.lat, names(stations))
  }
  if (class(obs.staid.time) == "character") {
    obs.staid.time <- match(obs.staid.time, names(obs))
  }
  ids <- sort(unique(stations[,stations.staid.lon.lat[1]]))
  
  if (class(obs[, obs.staid.time[1]]) != class(stations[, stations.staid.lon.lat[1]])) {
    stop("Classes of ID column in obs and stations is different!")
  }
  
  time <- unique(obs[, obs.staid.time[2]])
  time <- sort(time)# as.POSIXlt(sort(time))
  
  nt <- length(time) # num of dates
  ns <- length(ids) # num of stations
  
  
  tempdf <- data.frame(rep(ids,nt), rep(time,each=ns)) 
  names(tempdf) <- names(obs)[obs.staid.time]
  
  #   require(plyr)
  data <- plyr::join(tempdf, obs, type = "left", match = "first")
  
  # sort like 1st station 1st date, 2nd stations. 1st date ... # no need, tempdf is already sorted
  # data <- data[ order( data[, 1]), ] # sort by station id, it is always 1!
  # data <- data[ order( data[, 2]), ] # sort by time, it is always 2!
  row.names(data) <- 1:length(data[,1])
  
  # system.time( merge(tempdf,obs, all=TRUE) )
  # system.time(join(tempdf,obs) )
  # join is 2 x faster
  ids <- data.frame(staid=ids)
  # ids <- ids[ order( ids[, 1]), ] # already done with ids
  # ids <- as.data.frame(ids)
  names(ids) <- names(stations) [ stations.staid.lon.lat[1] ]
  
  st <- join(ids, stations)
  names(st)[ stations.staid.lon.lat[2:3] ] <- c('lon', 'lat')
  coordinates(st) <-~ lon +lat
  st@proj4string <- crs
  
  data <- as.data.frame(data[,-c(1,2)])
  names(data)= names(obs)[-obs.staid.time]
  
  if (is.null(delta) && length(time)==1){
    endTime <- as.POSIXct(time + 86400)
    stfdf <-STFDF(st, time , data, endTime)
  } else if (is.null(delta) && length(time)!=1) {
    stfdf <-STFDF(st, time , data)
  } else {
    endTime <- as.POSIXct(time + delta)
    stfdf <-STFDF(st, time , data, endTime)
  }

  # count NAs per stations
  
  bools2 <- c()
  for (i in 1:ncol(stfdf@data)){
    bools2 <- cbind(bools2, apply(matrix(stfdf@data[, i],
                           nrow=length(stfdf@sp),byrow=F), MARGIN=1,
                    FUN=function(x) sum(is.na(x))))
  }
  bools2 <- apply(bools2, 1, sum)
  # remove all NA
  stfdf=stfdf[bools2!=(nt*ncol(stfdf@data)), ,drop=F]
  
  row.names(stfdf@sp) <- 1:nrow(stfdf@sp)
  
  return(stfdf)
  
}
