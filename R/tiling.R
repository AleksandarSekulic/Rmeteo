tiling <- function(rast="", # path to grid file in SAGA format / raster formats 
                   tilesize=500, # in cells nx= ny
                   overlapping=50, # in cells
                   aspoints= NA, # "sf" # "sp"
                   asfiles=FALSE,                   
                   tilename="tile", # prexif to be given to tile names
                   tiles_folder='tiles', # resulting folder
                   parallel.processing=FALSE,
                   cpus=6,
                   ...) # filetype="GTiff", overwrite = T, datatype = "INT2S"
{
  
  if(inherits(rast, "SpatialPixelsDataFrame")) {
    r=rast(rast)
    ext(r) <- c(bbox(rast)[1,1], bbox(rast)[1,2], bbox(rast)[2,1], bbox(rast)[2,2])
  } else if (inherits(rast, "SpatRaster")){
    # do nothing
    r <- rast
  } else {r= rast(rast)}
  
  if (length(tilesize) == 1) {
    tilesize <- rep(tilesize, 2)
  }
  if (length(overlapping) == 1) {
    overlapping <- rep(overlapping, 2)
  }
  
  bb <- ext(r)
  bb <- cbind(c(bb$xmin,bb$xmax),c(bb$ymin,bb$ymax))
  
  cel1 <- res(r) [1]
  cel2 <- res(r) [2]
  
  bb[1,]=bb[1,] - overlapping[1]*cel1
  bb[2,]=bb[2,] + overlapping[2]*cel2
  
  step_x=tilesize[1]*cel1 + overlapping[1]*cel1
  step_y=tilesize[2]*cel2 + overlapping[2]*cel2
  nlon=ceiling(diff(bb[,1]) / step_x )
  nlat=ceiling(diff(bb[,2]) / step_y )
  
  p.l <- expand.grid(KEEP.OUT.ATTRS=FALSE,
                     lonmin=seq(bb[1,1],bb[1,1]+(nlon-1)*step_x, by=step_x),
                     latmin=seq(bb[1,2],bb[1,2]+(nlat-1)*step_y, by=step_y)) # - overlapping[2]*cel2    
  p.l$lonmin[p.l$lonmin<ext(r)$xmin] = ext(r)$xmin
  p.l$latmin[p.l$latmin<ext(r)$ymin] = ext(r)$ymin
  p.u <- expand.grid(KEEP.OUT.ATTRS=FALSE,
                     lonmax=seq(bb[1,1]+step_x, bb[1,1]+nlon*step_x, by=step_x),
                     latmax=seq(bb[1,2]+step_y   ,bb[1,2]+nlat*step_y, by=step_y)) + overlapping[2]*cel2  
  p.u$lonmax[p.u$lonmax>ext(r)$xmax] = ext(r)$xmax
  p.u$latmax[p.u$latmax> ext(r)$ymax] = ext(r)$ymax
  
  ptiles <- cbind(p.l, p.u)
  
  poligoni=as.list(rep(NA,length(ptiles$lonmin)))
  
  for(i in 1:length(ptiles$lonmin)) {
    x <- rbind(ptiles$lonmin[i], ptiles$lonmax[i],ptiles$lonmax[i],ptiles$lonmin[i],ptiles$lonmin[i])
    y <- rbind(ptiles$latmin[i], ptiles$latmin[i], ptiles$latmax[i], ptiles$latmax[i],ptiles$latmin[i])
    poligoni[[i]] <- sf::st_polygon(list(cbind(x,y)))
    # poligoni[[i]] <- Polygons( list(Polygon(cbind(x,y))) ,i)  
  }
  
  poll <- sf::st_multipolygon(poligoni)
  # poll <- SpatialPolygons(poligoni)
  # plot(poll)
  
  tiles=as.list(rep(NA,length(poll)))
  
  if(parallel.processing){
    if(!sfParallel()){
      sfInit(parallel = TRUE, cpus =cpus) 
      sfst = FALSE
    } else {sfst=TRUE}
    
    sfLibrary(package="terra", character.only=TRUE)
    sfExport( "poll" )
    sfExport( "r" )
    sfExport( "tiles_folder" )
    sfExport( "tilename" )
    # sfExport( "format" )
    if(asfiles){
      if(!dir.exists(tiles_folder)){
        dir.create(tiles_folder)
      }
      filename <- paste(tiles_folder,"/",tilename,i,sep="")
    } else {
      filename <- ""
    }
    tiles <- sfLapply(1:length(poll), function(i)  {
      terra::crop(r, poll[[i]], filename = filename, ...)
      # writeRaster(tiles[[1]], filename= paste(tiles_folder,"/",tilename,i,".",format,sep=""), ...) } )
    })
  } else {
    if(asfiles){
      if(!dir.exists(tiles_folder)){
          dir.create(tiles_folder)
      }
      filename <- paste(tiles_folder,"/",tilename,1:length(poll),sep="")
    } else {
      filename <- rep("", length(poll))
    }
    for(i in 1:length(poll)){
      tiles[i] <- terra::crop(r, poll[[i]], filename = filename[i], ...) # filetype=filetype, overwrite=overwrite)
      # writeRaster(tiles[i][[1]], filename= paste(tiles_folder,"/",tilename,i,".",format,sep=""), ...)
    }
  }
  
  if(parallel.processing) { sfStop() }
  
  if(!is.na(aspoints)) {
    rem <- sapply(tiles, function(i) all(is.na(values(i))))# i@data@values) ) )
    tiles = tiles[!rem]
    # tiles= lapply(tiles, function(i) rasterToPoints(i,spatial=TRUE) )
    tiles = lapply(tiles, function(i) as.points(i))
    if (aspoints == "terra") {
      # do nothing
    } else if (aspoints == "sf") {
      tiles = lapply(tiles, function(i) st_as_sf(i))
    } else if (aspoints == "sp") {
      tiles = lapply(tiles, function(i) as(i, "Spatial"))
    }
  }
  return(tiles)
}

  

  

