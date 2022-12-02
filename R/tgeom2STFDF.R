temp_geom<-function(day,
                    fi,
                    variable="mean",
                    ab=NULL){
  
  if(is.null(ab)){
    if(variable=='min'){
      a <- 24.16453  
      b <-  -15.71751
    } else if(variable=='max'){
      a <- 37.03043  
      b <- -15.43029
    }else{ # mean
      a <- 30.419375
      b <- -15.539232
      if(variable!='mean'){
        warning("variable argument is not valid! 'mean' is used.")
      }
    }
  } else{
    a <- ab[1]
    b <- ab[2]
  }
  
  f=ifelse(fi==0,1e-10,fi)
  costeta= cos( (day-18 )*pi/182.5 +2^(1-sign(fi) ) *pi) 
  cosfi = cos(fi*pi/180 )
  A=cosfi
  B= (1-costeta ) * abs(sin(fi*pi/180 ) )
  x=a*A + b*B 
  return(x)}



tgeom2STFDF <- function(grid,
                        time,
                        variable='mean',
                        ab=NULL) {
  
  if(class(grid) == 'SpatialGrid' | class(grid) == 'SpatialGridDataFrame') {
  grid <- as(grid,'SpatialPixels') }
  
  if(!is.na(grid@proj4string@projargs )){ grid1 <- spTransform(grid, CRS('+proj=longlat +datum=WGS84') ) } else { grid1 <- grid}
  
  time <- as.POSIXlt(sort(time))
  day<- as.numeric( strftime(time, format = "%j") )
  
  tg<-lapply(day, function(i) temp_geom(i,grid1@coords[,2],ab) )
  tg=do.call('cbind',tg)
  tg=as.vector(tg)
  tg=data.frame(temp_geo=tg)
  
  endTime=time+86400

  res <- STFDF(grid,time,data=tg,endTime=endTime)
  
  return(res)
}
