library(meteo)
library(sp)
library(sf)
library(spacetime)
library(sftime)
library(terra)
library(raster)
library(gstat)
library(plyr)
library(xts)
library(snowfall)
library(doParallel)
library(CAST)


#################### tgeom2STFDF ####################
# sf
demo(meuse, echo=FALSE)
meuse <- meuse[complete.cases(meuse@data),]
locations = st_as_sf(meuse, coords = c("x", "y"), crs = 28992, agr = "constant") # sf

grid = locations
time = seq(as.POSIXct("2011-07-01"), as.POSIXct("2011-07-02"), by="day")
variable = "mean"
ab=NULL

# sftime
data(dtempc_ogimet)
dtempc <- dtempc[complete.cases(dtempc),]
# summary(dtempc)
dtempc.sf <- st_as_sf(dtempc, coords = c("lon", "lat"), crs = 4326, agr = "constant")
# dtempc.sf <- st_transform(dtempc.sf, 32634)
dtempc.sftime <- st_sftime(dtempc.sf)
class(dtempc.sftime)
st_crs(dtempc.sftime)
dtempc.sftime = dtempc.sftime[dtempc.sftime$time %in% unique(dtempc.sftime$time)[1:2], ]
locations = st_drop_time(dtempc.sftime)

grid <- locations[!duplicated(locations$geometry), ]
time = unique(st_time(dtempc.sftime))
variable = "mean"
ab=NULL

# SpaVector
grid = terra::vect(meuse)
time = seq(as.POSIXct("2011-07-01"), as.POSIXct("2011-07-02"), by="day")
variable = "mean"
ab=NULL

# SpatRaster
grid <- terra::rast(meuse.grid)
class(grid)
time = seq(as.POSIXct("2011-07-01"), as.POSIXct("2011-07-02"), by="day")
variable = "mean"
ab=NULL

temp_geo <- tgeom2STFDF(grid = grid,
                        time = time,
                        variable=variable)

#################### tilling ####################

demo(meuse, echo=FALSE)
rast <- terra::rast(meuse.grid[, "dist"])
rast <- meuse.grid[, "dist"]
# 104 xx 78 x 5
# rast="" # path to grid file in SAGA format / raster formats 
tilesize=20 # in cells nx= ny
overlapping=5 # in cells
aspoints= "sf" # "sp", "terra" # NA
asfiles=TRUE # FALSE                   
tilename="tile" # prexif to be given to tile names
tiles_folder=paste(getwd(),'tiles',sep='/') # resulting folder
parallel.processing=FALSE
cpus=6
filetype="GTiff"
overwrite = T
datatype = "INT2S"

tiles <- tiling(rast=rast, # path to grid file in SAGA format / raster formats 
                tilesize=tilesize, # in cells nx= ny
                overlapping=overlapping, # in cells
                aspoints= aspoints, # "sf" # "sp"
                asfiles=asfiles,                   
                tilename=tilename, # prexif to be given to tile names
                tiles_folder=tiles_folder, # resulting folder
                parallel.processing=parallel.processing,
                cpus=cpus,
                filetype=filetype,
                overwrite = overwrite)
                # datatype = datatype)
plot(tiles[[1]])

#################### meteo2STFDF ####################

data(dtempc_ogimet) 
data(stations_ogimet)
# stations <- unique(dtempc[, c(1,4,3,5,2,8,9)])
lonmin=18 ;lonmax=22.5 ; latmin=40 ;latmax=46
serbia = point.in.polygon(stations$lon, stations$lat, c(lonmin,lonmax,lonmax,lonmin), 
                          c(latmin,latmin,latmax,latmax))
st= stations[ serbia!=0, ] # stations in Serbia approx.
obs = dtempc[, c(1,6,7,10,11,12)]
stations = st
obs.staid.time=c("staid", "time")
stations.staid.lon.lat=c(1,2,3)
crs= CRS('+proj=longlat +datum=WGS84')
delta=NULL

# create STFDF
temp <- meteo2STFDF(obs = obs,
                    stations = stations,
                    crs = crs)

#################### rm.dupl ####################

data(dtempc_ogimet) 
data(stations_ogimet)
# stations <- unique(dtempc[, c(1,4,3,5,2,8,9)])

# create STFDF
temp <- meteo2STFDF(dtempc[, c(1,6,7,10,11,12)],
                    stations,
                    crs = CRS('+proj=longlat +datum=WGS84')) # create STFDF from 2 data frames
nrow(temp@sp) # number of stations before removing dupl.
temp2 <-rm.dupl(temp, zcol = 1, zero.tol = 50)
nrow(temp2@sp) # number of stations after

#################### rfilltimegaps / rfillspgaps ####################
# sp gaps
data(nlmodis20110704)
data(nlmodis20110712)
data(NLpol)

# terra
nlmodis20110704 <- terra::rast(nlmodis20110704)
nlmodis20110712 <- terra::rast(nlmodis20110712)

# SpaVector
NLpol = vect(NLpol)
crs(NLpol) <- "epsg:4326"
# # sf
# NLpol <- st_as_sf(NLpol) #, crs = st_crs(4326))

plot(nlmodis20110712)
plot(NLpol, add=T)

# fill spatial gaps
nlmodis20110704 <- rfillspgaps(nlmodis20110704, NLpol)
nlmodis20110712 <- rfillspgaps(nlmodis20110712, NLpol)
plot(nlmodis20110712)
summary(values(nlmodis20110704))
summary(values(nlmodis20110712))

# time
nlmodis20110704 <- as(raster(nlmodis20110704),"SpatialPixelsDataFrame")
names(nlmodis20110704)='m1'
nlmodis20110712 <- as(raster(nlmodis20110712),"SpatialPixelsDataFrame")
names(nlmodis20110712)='m2'

nlmodis20110704@data <- cbind(nlmodis20110704@data, nlmodis20110712@data)

df<-reshape(nlmodis20110704@data , varying=list(1:2), v.names="modis",direction="long", 
            times=as.Date(c('2011-07-04','2011-07-12')), ids=1:dim(nlmodis20110704)[1])

stMODIS<- STFDF(as( nlmodis20110704, "SpatialPixels"), 
                time= as.Date(c('2011-07-04','2011-07-12')), 
                data.frame(modis=df[,'modis']))
stMODIS <- rfilltimegaps(stMODIS)
# stplot(stMODIS, col.regions=bpy.colors())

