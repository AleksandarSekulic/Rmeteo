\name{tiling}
\alias{tiling}
\title{Tiling \link[raster]{raster} or \link[sp]{Spatial-class} Grid or Pixels object}
\description{Tiling \link[raster]{raster} or \link[sp]{Spatial-class} Grid or Pixels (data frame) object to smaller parts with optional overlap.
}
\usage{
tiling(filename,tilesize=500,overlapping=50, aspoints= FALSE, 
       asfiles=FALSE,tilename="tile", format="GTiff", 
       tiles_folder=paste(getwd(),'tiles',sep='/'), parallel.processing=FALSE, cpus=6)
}
\arguments{
  \item{filename}{raster object, SpatialPixels* object, SpatialGrid* object or file path of raster object stored on the disk (can be read via \link[rgdal]{readGDAL}), for more details see \link[raster]{raster}. The  resolution of object should be the same in x and y direction.
}
  \item{tilesize}{tile size in cells in x direction. nx = ny is assumed, total number of tile cells is \code{tilesize} by \code{tilesize}. }
  \item{overlapping}{overlapping cells in each direction}
  \item{aspoints}{if TRUE tiles are in form of \code{SpatialPointsDataFrame} }
  \item{asfiles}{if TRUE tiles are stored on local drive as raster objects}
  \item{tilename}{prefix given to file names}
  \item{format}{file format, see \link[raster]{writeRaster} }
  \item{tiles_folder}{folder to be created for tiles storage}
  \item{parallel.processing}{if TRUE parralel processing is performed via  \link[snowfall]{sfLapply}
}
  \item{cpus}{number of proccesing units}
}
\value{ The list of tiles in raster format or in \code{SpatialPointsDataFrame} format if \code{aspoints=TRUE}
}
\author{Milan Kilibarda \email{kili@grf.bg.ac.rs}  }
\seealso{
\code{\link{pred.strk}}
}
\examples{
data(regdata) 
str(regdata@sp) # DEM and TWI data for Serbia at 1 km resolution
dem=regdata@sp['dem']
library(sp)
spplot(dem,  col.regions=bpy.colors() )

str(dem)
# tiling dem in tiles 250x250 with 25 cells overlap
tiles = tiling(dem,tilesize=250,overlapping=25,aspoints=TRUE)
# number of tiles
length(tiles)

image(dem)
points(tiles[[1]] , pch='-' ,col ='green')
points(tiles[[2]], pch='.')

str(tiles[[1]])

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

